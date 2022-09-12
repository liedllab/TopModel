"""Pymol script to visualize the structure checked by app.py"""
import pickle
from pymol import cmd # pylint: disable=import-error

def main_pymol():
    """Startup script to work with app.py"""
    with open('.pdbear.temp', 'rb') as temp:
        path = pickle.load(temp)
        colors = pickle.load(temp)
        data = pickle.load(temp)

    cmd.load(path)
    cmd.color('white')
    cmd.set('cartoon_transparency', 0.3)
    cmd.hide('lines')
    cmd.show('cartoon')
    for k, values in data.items():
        if values:
            if not k in ('D', 'L'):
                sel_str = " or ".join(["resid " + str(n) +  "-" + str(n+1) for n in values])
            else:
                sel_str = " or ".join(["resid " + str(n) for n in values])

            try:
                cmd.color(colors[k], sel_str)
                cmd.show("sticks", sel_str)
            except KeyError:
                pass
            cmd.select(k, sel_str)

    cmd.color('atomic', 'not elem C')
    cmd.deselect()


if __name__ == 'pymol':
    main_pymol()
