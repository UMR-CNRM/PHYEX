# -*- coding: utf-8 -*-
"""
Created on Wed Mar 30 22:32:56 2022

@author: Quentin
"""

def count_blank(text):
    count=0
    for t in text:
        if t == " ":
            count=count+1
        else:
            break
    return count
        
def check_indent(n,indent_score,worktext):
    nblank=count_blank(worktext)
    rawtext=worktext[nblank:]
    correcttext=' '*((indent_score)*1) + rawtext
    return correcttext

filein = open('mode_tke_eps_sources.f90','r')
fileout  = open('new.F90','w')
contentbyline = filein.readlines()

ncurrline=0
indent_score=0
for i in contentbyline:
    #Si la ligne est un commentaire (commence par !) : next (ne pas indenter)
    #
    # ligne de calcul sans ajout d'indentation (DO, ENDDDO, IF, ENDIF)
    if indent_score > 0 and ("DO J" not in i and "ENDDO" not in i):
        textwrite=check_indent(ncurrline,indent_score,i)
    elif "DO J" in i:
        textwrite=check_indent(ncurrline,indent_score,i)
        indent_score = indent_score+1
    elif "ENDDO" in i:
        indent_score = indent_score-1
        textwrite=check_indent(ncurrline,indent_score,i)
    #Ajouter le case IF et ENDIF qui ajoute de lindentation
    else:
        textwrite=i
    ncurrline=ncurrline+1
    fileout.writelines(textwrite)

fileout.close()
filein.close()
    