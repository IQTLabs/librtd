# Package

version       = "0.1.0"
author        = "Benjamin Lee"
description   = "Generalized k-mer return time distribution calculation"
license       = "Apache-2.0"
srcDir        = "src"
installExt    = @["nim"]
bin           = @["librtd"]



# Dependencies

requires "nim >= 1.0.4"
requires "docopt >= 0.6.8"
requires "progress >= 1.1.1"