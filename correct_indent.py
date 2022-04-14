#!/usr/bin/env python3
def detectMNH_expand(f):
    fin = open(f,'r')
    fout  = open(f+'_EXPAND','w')
    contentbyline = fin.readlines()
    
    for i in contentbyline:
        if "!$mnh_expand" in i:
            fout.writelines("! $MNH EXPAND$ !\n")
            fout.writelines(i)    
        elif "!$mnh_end_expand" in i:
            fout.writelines(i)
            fout.writelines("! $MNH END EXPAND$ !\n")
        else: 
            fout.writelines(i)
    
    fout.close()
    fin.close()

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

def firstcharacnotblank(string):
    for i in string:
        if i == " ": pass
        else: return i

def correct_indent(f):
    fin = open(f,'r')
    fout  = open(f+'_CORRECT_INDENT','w')
    contentbyline = fin.readlines()
    
    ncurrline=0
    indent_score=0
    expand_score=0
    for i in contentbyline:
        if "! $MNH EXPAND$ !" in i:
            expand_score = expand_score + 1
            textwrite = ""
        elif "! $MNH END EXPAND$ !" in i:
            expand_score = expand_score - 1
            textwrite = ""
        # Correct indentation only in between $MNH EXPAND$ and $MNH END EXPAND$
        elif expand_score >= 1:
            # Do not indent comment lines within mnh_expand
            if firstcharacnotblank(i) == "!":
                textwrite=i
            elif indent_score > 0 and ("DO J" not in i and "ENDDO" not in i and "END IF" not in i and "ENDIF" not in i and "THEN" not in i and "ELSE" not in i):
                textwrite=check_indent(ncurrline,indent_score,i)
            elif "ELSE" in i:
                indent_score = indent_score+1        
                textwrite=check_indent(ncurrline,indent_score,i)
                indent_score = indent_score-1
            elif "THEN" in i:
                textwrite=check_indent(ncurrline,indent_score,i)
                indent_score = indent_score+1        
            elif "DO J" in i:
                textwrite=check_indent(ncurrline,indent_score,i)
                indent_score = indent_score+1
            elif ("ENDDO" in i) or ("END DO" in i):
                indent_score = indent_score-1
                textwrite=check_indent(ncurrline,indent_score,i)
            elif ("END IF" in i) or ("ENDIF" in i):
                indent_score = indent_score-1
                textwrite=check_indent(ncurrline,indent_score,i)       
            else:
                textwrite=i
        else: # not EXPAND lines nor within mnh_expand
            if "THEN" in i:
                indent_score = indent_score+1        
            elif "DO J" in i:
                indent_score = indent_score+1
            elif ("ENDDO" in i) or ("END DO" in i):
                indent_score = indent_score-1
            elif ("END IF" in i) or ("ENDIF" in i):
                indent_score = indent_score-1
            textwrite=i
        ncurrline=ncurrline+1
        fout.writelines(textwrite)
    
    fout.close()
    fin.close()

if __name__ == "__main__":
    import argparse
    import sys
    parser = argparse.ArgumentParser(description='Detecte les balises !$mnh_expand et !$mnh_end_expand et ajoute un commentaire avant et après pour corriger l\'indentation par la suite avec correct_indentation.py' )
    value = argparse.ArgumentParser()
    parser.add_argument('file1', metavar='file1', type=str, help="file to check")
    parser.add_argument('action', metavar='action', type=str, help="action (indent or detect")   
    args = parser.parse_args()
    if "indent" in args.action:
        sys.exit(correct_indent(args.file1))
    elif "detect" in args.action:
        sys.exit(detectMNH_expand(args.file1))
    else:
        sys.exit("Error : action should be indent or detect. Nothing has been done")
