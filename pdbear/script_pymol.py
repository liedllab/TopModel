import pickle
from pymol import cmd
import sys
import os

def main_pymol():
    with open('.temp', 'rb') as f:
        path = pickle.load(f)
        colors = pickle.load(f)
        data = pickle.load(f)
    
    cmd.load(path)
    cmd.color('white')
    cmd.set('cartoon_transparency', 0.4)
    cmd.hide('lines')
    cmd.show('cartoon')

    for k, values in data.items():
        if values:
            if not k=='D':
                sel_str = " or ".join(["resid " + str(n) +  "-" + str(n+1) for n in values])
            else:
                sel_str = " or ".join(["resid " + str(n) for n in values])
            print(sel_str)
            print(k)

            cmd.color(colors[k], sel_str)
            cmd.show("sticks", sel_str)
            cmd.select(k, sel_str)
    
    cmd.color('atomic', 'not elem C')
    os.remove('.temp')


if __name__ == 'pymol':
    main_pymol()
