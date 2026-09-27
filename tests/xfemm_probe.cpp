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
int main(int argc, char **argv) {
    if (argc != 4) return 2;
    FPProc reader;
    if (!reader.OpenDocument(argv[1])) return 3;
    reader.Smooth = false;
    CMPointVals p;
    if (!reader.GetPointValues(std::strtod(argv[2], nullptr), std::strtod(argv[3], nullptr), p)) return 4;
    std::cout << std::setprecision(17) << "VALUES " << p.A.re << " " << p.B1.re << " " << p.B2.re
              << " " << p.H1.re << " " << p.H2.re << "\n";
    return 0;
}
