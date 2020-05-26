# Package

version       = "0.1.0"
author        = "Benjamin Lee"
description   = "Generalized k-mer return time distribution calculation"
license       = "Apache-2.0"
srcDir        = "librtd"
installExt    = @["nim"]
skipFiles     = @["cli.nim"]


# Dependencies

requires "nim >= 1.0.4"
requires "progress >= 1.1.1"
requires "nimpy >= 0.1.0"