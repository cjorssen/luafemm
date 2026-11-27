#!/usr/bin/env python3
# SPDX-License-Identifier: LPPL-1.3c
# Copyright (C) 2026 Christophe Jorssen
# Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
# LPPL maintenance status: maintained
# This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
"""Compile the manual's configuration fragments in explicit physical contexts.

Rendered figures execute inside the manual. These checks supply the surrounding
picture/axis or previously captured model required by reference fragments; the
fragment bytes themselves are extracted unchanged from the displayed source.
A failed calculation fails the documentation build. No output from this checker
is imported into the manual.
"""
from pathlib import Path
import re
import json
import hashlib


def prepare(root, output):
    snippets={}
    for p in sorted((Path(__file__).resolve().parent/'sections').glob('*.tex')):
        for n,m in enumerate(re.finditer(r'\\begin\{codeexample\}\[code only\](.*?)\\end\{codeexample\}',p.read_text(),re.S),1):
            snippets[f'{p.stem}:{n}']=m[1].strip()
    def s(name): return snippets[name]
    def picture(code, options='femm/problem={xmin=-110,xmax=110,ymin=-95,ymax=115,mesh size=8,depth=20}'):
        return '\\begin{tikzpicture}['+options+']\n'+code+'\n\\end{tikzpicture}\n'
    def axis(code):
        return '\\begin{tikzpicture}\\begin{axis}[width=10cm,height=5cm,ymin=-1,ymax=1]\n'+code+'\n\\end{axis}\\end{tikzpicture}\n'
    styles='\n'.join(s('tutorial:'+str(i)) for i in [2,4,6,12])
    coil=s('tutorial:5')
    mag=r'\filldraw[core] (-40,-30) rectangle (40,-15);'+s('tutorial:13')
    circuit=r'\path[femm/region={material=copper,circuit=drive,turns=10,winding fill=.6}] (-12,-5) rectangle (-8,5);\path[femm/region={material=copper,circuit=drive,turns=-10,winding fill=.6}] (8,-5) rectangle (12,5);'
    bias=picture(circuit,'femm/problem={xmin=-30,xmax=30,ymin=-25,ymax=25,mesh size=5,depth=20},femm/circuit={name=drive,current=5},femm/operating point={name=bias}')
    biasrod=picture(circuit,'femm/problem={xmin=-30,xmax=30,ymin=-25,ymax=25,mesh size=5,depth=20},femm/circuit={name=drive,current=5},femm/operating point={name=rod-low}')
    slabs=s('harmonic:1')
    slabbody=s('harmonic:2').split('% In a PGFPlots axis:')[0]
    ndtstyle=s('tutorial-ndt:1')
    ndtbody=s('tutorial-ndt:2')
    # The documented force/capture fragments require complete physical models.
    # These independent harnesses use the corresponding downloadable geometry;
    # the checker discards their output, never importing it into the manual.
    sliding=r'\input{build/doc/bodies/sliding-band.tex}'
    force=r'\input{build/doc/bodies/forces.tex}'
    contexts={}
    for name,code in snippets.items():
        section,n=name.split(':');n=int(n)
        if code.startswith('local '):
            # Preserve displayed bytes; TeX joins directlua lines, so only the
            # execution copy drops Lua line comments. Never touch TikZ paths.
            code=re.sub(r'--[^\n]*','',code)
        rendered=None
        if section in ('getting-started','tutorial') and (section=='getting-started' or n==1):
            # These are loading skeletons. Native format versions are checked
            # separately below, with a small page so fresh output is required.
            contexts[name]={'format':'latex' if '\\documentclass' in code else 'context' if '\\starttext' in code else 'plain','source':code}
            continue
        if section=='tutorial':
            if n in [2,3,4,6,12,17]:rendered=code
            elif n in [5,7,8,9,10,13,14,15]:rendered=styles+picture((coil if n in [8,9] else '')+code)
            elif n==11:rendered=styles+picture(coil+s('tutorial:10'))+code
            elif n==16:rendered=styles+picture(mag+s('tutorial:15'))+code
            elif n==18:
                # A two-part key list, explicitly separated in the displayed text.
                parts=code.split('% In the style definition:')
                rendered='\\tikzset{'+parts[0]+'}\n\\tikzset{'+parts[1]+'}'
            elif n==19:
                rendered=styles+picture(coil+r'\path[femm/profile={name=coil-coarse}] (-48,31.5)--(48,31.5);\path[femm/profile={name=coil-fine}] (-48,31)--(48,31);')+axis(code)
        elif section=='axisymmetric':
            if n==1:rendered=code
            elif n==2:rendered=s('axisymmetric:1')+picture(code,'solenoid')
            elif n==3:
                a,b=code.split('% In an ordinary PGFPlots axis:')
                rendered=s('axisymmetric:1')+picture(s('axisymmetric:2')+a,'solenoid')+axis(b)
            elif n==4:rendered=s('axisymmetric:1')+picture(s('axisymmetric:2')+r'\path[femm/profile={name=name}] (0,-10)--(0,10);'+code,'solenoid')
        elif section=='operating':
            if n in [1,5]:rendered=code
            elif n in [2,3]:rendered=bias+code
            elif n==4:
                measures=''.join(r'\femmmeasure[quantity=circuit,circuit=drive,name='+v+',parameter='+str(i)+']' for i,v in enumerate(['small','small-a','small-b','small-c']))
                rendered=bias+picture(measures,'femm/problem={operating point=bias,permeability=incremental,circuit currents={drive=1}}')+axis(code)
            elif n==6:rendered=bias+'\\directlua{local model=luafemm.current\n'+code+'\n}'
        elif section=='biased-harmonic':
            if n==1:rendered=biasrod+code
            elif n==2:rendered=r'\tikzset{test/.style={femm/problem={'+code+'}}}'
            elif n==3:
                measures=''.join(r'\femmmeasure[quantity=circuit,circuit=drive,name='+v+',parameter='+str(i)+']' for i,v in enumerate(['terminal','port-low-5','port-low-10','port-low-20']))
                rendered=biasrod+picture(measures,'femm/problem={operating point=rod-low,permeability=incremental,frequency=200,circuit currents={drive=1}}')+axis(code)
            elif n==4:rendered=biasrod+'\\directlua{local bias=require("luafemm-operating").capture(luafemm.current)\n'+code+'\n}'
        elif section=='sliding':
            if n==1:rendered=sliding+code
            elif n==2:rendered=sliding+picture(code,'femm/problem={sliding band=machine,rotor angle=60}')
            elif n==3:
                # Reuse the uncaptured declaration body, without its capture key.
                body=(output/'bodies/sliding-band.tex').read_text()
                m=re.search(r'\\begin\{tikzpicture\}\[femm/problem=.*?\\end\{tikzpicture\}',body,re.S)
                if not m:raise RuntimeError('Missing sliding declaration')
                declaration=re.sub(r',\s*femm/sliding band=\{[^}]+\}', '',m[0])
                rendered=body.split(r'\centering')[0]+declaration+'\\directlua{local model=luafemm.current\n'+code+'\n}'
        elif section=='tutorial-ndt':
            if n==1:rendered=code
            elif n in [2,3,5]:rendered=ndtstyle+picture(('' if n==2 else ndtbody)+code,'ndt/model')
            elif n==4:
                a,b=code.split('% Later, inside a PGFPlots axis:');rendered=ndtstyle+picture(ndtbody+a,'ndt/model')+axis(b)
            elif n==6:rendered=ndtstyle+r'\def\Position{0}\def\Case{test}\def\Index{0}'+picture(code,'ndt/model')
            elif n==7:rendered=ndtstyle+picture(ndtbody+r'\femmmeasure[name=deep-10,quantity=circuit,circuit=probe]\femmmeasure[name=sound-10,quantity=circuit,circuit=probe]','ndt/model')+code
        elif section=='harmonic':
            if n==1:rendered=code
            elif n==2:
                a,b=code.split('% In a PGFPlots axis:');rendered=slabs+picture(a,'conducting slab')+axis(b)
            elif n==3:
                a,b=code.split('% In this picture:');rendered=a+picture(b,'induction')
            elif n==4:rendered=slabs+picture(slabbody+code,'conducting slab')
            elif n==5:rendered=slabs+r'\tikzset{femm/materials/solid copper/.style={library=copper}}'+picture(slabbody.replace('slab copper','solid copper')+code,'conducting slab')
        elif section=='currents':
            a,b=s('currents:1').split('% In a picture with this style:')
            if n==1:rendered=a+picture(b,'induction experiment')
            elif n==2:
                c,d=code.split('% Later, in a PGFPlots axis:');rendered=a+picture(b+c,'induction experiment')+axis(d)
            elif n==3:
                rendered=picture(r'\path[femm/region={material=copper}] (1,-2) rectangle (3,2);'+code,'femm/problem={type=axisymmetric,frequency=50,xmax=8,ymin=-5,ymax=5,mesh size=1}')
            elif n==4:rendered=a+picture(b+code,'induction experiment')
        elif section=='tutorial-machines':
            if n==1:rendered=code
            elif n==2:rendered=s('tutorial-machines:1')+picture(code)
            elif n==3:rendered=s('tutorial-machines:1')+picture(s('tutorial-machines:2'))+'\\begin{tikzpicture}'+code+'\\end{tikzpicture}'
            elif n in [4,5]:rendered=r'\tikzset{variant/.style={femm/stator={'+(s('tutorial-machines:4')+',' if n==5 else '')+code+'}}}'+picture(r'\node[variant] {};')
        elif section=='postprocessing':
            rendered=picture(r'\node[femm/u electromagnet={gap=3,gap mesh size=1,ampere turns=2000,core material=armco-pure-iron,armature material=armco-pure-iron},anchor=core center] {};'+(r'\path[femm/selection={name=armature,material=armco-pure-iron}] (-42,31) rectangle (42,55);' if n==2 else '')+code)
        elif section=='ac-circuits':
            rendered=picture(r'\path[femm/region={material=copper,circuit=supply,winding=lower,turns=10,winding fill=.6}] (-4,-4) rectangle (-2,4);'+code,'femm/problem={frequency=50,depth=20,xmin=-10,xmax=10,ymin=-10,ymax=10,mesh size=2},femm/circuit={name=supply,current=1}')
        elif section=='components':rendered=picture(code)
        elif section=='boundaries':rendered=picture(code) if n==2 else code
        elif section in ['materials','losses','harmonic-material-models']:rendered=code
        elif section=='development':rendered='\\directlua{'+code+'}'
        if rendered is None:raise RuntimeError('No execution context for '+name)
        contexts[name]={'format':'fragment','source':rendered}
    checks=output/'snippets';checks.mkdir(parents=True,exist_ok=True)
    records={}
    for name,context in contexts.items():
        filename=name.replace(':','-')
        records[name]={'sha256':hashlib.sha256(snippets[name].encode()).hexdigest(),'format':context['format']}
        (checks/(filename+'.tex')).write_text(context['source']+'\n')
    (checks/'inventory.json').write_text(json.dumps(records,indent=2)+'\n')
    return contexts


def check(builder):
    """Execute every fragment and loading skeleton, plus the README example."""
    root, output = builder.ROOT, builder.BUILD / 'doc'
    contexts = prepare(root, output)
    preamble = r"""\documentclass{article}
\usepackage[margin=10mm]{geometry}
\usepackage{amsmath,tikz,pgfplots}
\usetikzlibrary{femm,calc,patterns}
\usepgfplotslibrary{femm}
\pgfplotsset{compat=1.18}
\tracinglostchars=3
\begin{document}
"""
    parts=[preamble]
    for name,context in contexts.items():
        if context['format'] != 'fragment':
            continue
        code=context['source'].replace('build/doc/bodies/', str(output/'bodies')+'/')
        parts.append('\\typeout{CHECK-SNIPPET: '+name+'}\n'
            + r'\begingroup\directlua{luafemm.tex.operating_points={};luafemm.tex.sliding_bands={}}'
            + r'\setbox0=\vbox{\hsize=180mm\hfuzz=1000pt\vfuzz=1000pt'+'\n'+code+'\n}\\endgroup\n')
    parts.append(r'Snippet execution checks.\end{document}')
    source=output/'snippet-checks.tex';source.write_text('\n'.join(parts))
    builder.run(['lualatex','--no-shell-escape','-interaction=nonstopmode','-halt-on-error',source.name],
                cwd=output,log=output/'snippet-checks-console.log')
    builder.check_log(output/'snippet-checks.log')
    log=(output/'snippet-checks.log').read_text()
    for name,context in contexts.items():
        if context['format']=='fragment':
            if 'CHECK-SNIPPET: '+name not in log:raise RuntimeError('Unexecuted fragment: '+name)
            continue
        folder=output/'snippets'/name.replace(':','-');folder.mkdir(parents=True,exist_ok=True)
        text=context['source'];fmt=context['format']
        ending={'latex':r'\end{document}','plain':r'\bye','context':r'\stoptext'}[fmt]
        text=text.replace(ending, 'Loading check.\n'+ending)
        (folder/'consumer.tex').write_text(text)
        (folder/'consumer.pdf').unlink(missing_ok=True)
        command=(builder.context_command() if fmt=='context' else
                 ['lualatex' if fmt=='latex' else 'luatex','--no-shell-escape','-interaction=nonstopmode','-halt-on-error'])
        builder.run(command+['consumer.tex'],cwd=folder,log=folder/'console.log')
        builder.check_log(folder/'consumer.log')
        if not (folder/'consumer.pdf').is_file():raise RuntimeError('Missing loading check PDF: '+name)
    readme=(root/'README.md').read_text()
    code=re.search(r'```tex\n(.*?)```',readme,re.S)[1]
    source=output/'readme-example.tex';source.write_text(code)
    builder.compile_tex(source,directory=output)
    (output/'snippets/passed.json').write_text(json.dumps({'count':len(contexts),'passed':True,'readme':hashlib.sha256((root/'README.md').read_bytes()).hexdigest()},indent=2)+'\n')
    print(f'PASS: all {len(contexts)} reference fragments/loading skeletons and the README example execute.')
