"""
Re-create the isentropic flow table as seen in
'Fundamentals of Aerodynamics' by J.D Anderson.
"""

import scipy as sp
from aerospacetoolbox import flowisentropic

def tablestr(v):
    e = sp.floor(sp.log10(v)) + 1
    m = v / 10**e
    s = ' $-$ ' if e < 0 else ' $+$ '
    return '%.4f' % m + s + '%02d' % sp.absolute(e)

M = sp.concatenate((sp.arange(0, 2, 0.02) + 0.02,
                    sp.arange(2, 5, 0.05) + 0.05,
                    sp.arange(5, 8, 0.1) + 0.1,
                    sp.arange(8, 20, 1) + 1,
                    sp.arange(20, 50, 2) + 2))

M, T, P, rho, area = flowisentropic(M=M)

with open('main.tex', 'w') as f:
    f.write('\\documentclass[a4paper,11pt]{article}\n')
    f.write('\\usepackage[a4paper, top=2cm]{geometry}\n')
    f.write('\\usepackage{longtable}\n')
    f.write('\\title{Appendix}\n')
    f.write('\\author{}\n')
    f.write('\\date{}\n')
    f.write('\\begin{document}\n')
    f.write('\\maketitle\n')
    f.write('\\begin{center}\n')
    f.write('\\section*{\\huge Isentropic Flow Properties}\n')
    f.write('\\setlength\\tabcolsep{10pt}\n')
    f.write('\\begin{longtable}{ c c c c c }\n')
    f.write('\\hline \\\\ $M$ & $\\frac{p_0}{p}$ & $\\frac{\\rho_0}{\\rho}$ & ' + \
            '$\\frac{T_0}{T}$ & $\\frac{A}{A^\\ast}$ \\\\ \\\\ \\hline \\\\ \n')
    f.write('\\endhead\n')
    f.write('\\\\ \\hline \\endlastfoot\n')
    for i in xrange(M.size):
        f.write(tablestr(M[i]) + ' & ' + \
                tablestr(1 / P[i]) + ' & ' + \
                tablestr(1 / rho[i]) + ' & ' + \
                tablestr(1 / T[i]) + ' & ' + \
                tablestr(area[i]) + ' \\\\\n')
        if (i+1) % 10 == 0 and not i == 0:
            f.write('\\\\\n')
    f.write('\\end{longtable}\n')
    f.write('\\end{center}\n')
    f.write('\\end{document}')

print "main.tex created and closed."
