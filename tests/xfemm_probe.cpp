// SPDX-License-Identifier: LPPL-1.3c
// Copyright (C) 2026 Christophe Jorssen
// Author and Current Maintainer: Christophe Jorssen <christophe.jorssen@gmail.com>
// LPPL maintenance status: maintained
// This file is part of luafemm. See LICENSE and LICENSES.md for its terms.
// Independent disk-reader oracle. This optional developer tool is not runtime code.
#include "fpproc.h"
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
int main(int argc, char **argv) {
    const bool post = argc == 5 && std::string(argv[1]) == "--post";
    const bool current = argc == 5 && std::string(argv[1]) == "--current";
    if (argc != 4 && !post && !current) return 2;
    const int offset = post || current ? 1 : 0;
    FPProc reader;
    if (!reader.OpenDocument(argv[1 + offset])) return 3;
    reader.Smooth = false;
    const double x = std::strtod(argv[2 + offset], nullptr);
    const double y = std::strtod(argv[3 + offset], nullptr);
    if (post) {
        const int element = reader.InTriangle(x, y);
        if (element < 0) return 4;
        for (auto &label : reader.blocklist) label.IsSelected = false;
        reader.blocklist[reader.meshelem[element].lbl].IsSelected = true;
        reader.WeightingScheme = 0;
        if (!reader.MakeMask()) return 5;
        std::cout << std::setprecision(17) << "POST";
        for (int type : {2, 17, 18, 19, 22})
            std::cout << " " << reader.BlockIntegral(type).re;
        std::cout << "\n";
        return 0;
    }
    CMPointVals p;
    if (!reader.GetPointValues(x, y, p)) return 4;
    std::cout << std::setprecision(17) << "VALUES " << p.A.re << " " << p.B1.re << " " << p.B2.re
              << " " << p.H1.re << " " << p.H2.re;
    if (current) std::cout << " " << p.Js.re * 1e6;
    std::cout << "\n";
    return 0;
}
