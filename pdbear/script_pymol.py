import sys

print(sys.argv)

#from pymol import cmd
#
#def pymol(path, data):
#    cmd.load(path)
#    cmd.color("white")
#    cmd.set("cartoon_transparency", 0.6)
#    cmd.select("D aminoacids", " or ".join([f"resid {n}" for n in data['D']]))
#    cmd.color("orange", "D aminoacids")
#    cmd.show("sticks", "D aminoacids")
#    
#    cmd.select("cis amides", " or ".join([f"resid {n}-{n+1}" for n in data['cis']]))
#    cmd.show("sticks", "cis amides")
#    cmd.color("red", "cis amides")
#    
#    cmd.select("to check", " or ".join([f"resid {n}-{n+1}" for n in data['?']]))
#    cmd.show("sticks", "to check")
#    cmd.color("yellow", "to check")
#
#cmd.extend("show", pymol)
#
#if __name__ == 'pymol':
#    pymol(*sys.argv[-2:])
