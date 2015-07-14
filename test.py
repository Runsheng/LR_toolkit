from utils import exe

print exe("ls -l")

from src.ssw_wrap import Aligner

Aligner.align("AAAAAAA","TTTTTTT")