#!/usr/bin/env python3

def detectMNH_expand(f):
# Adds MNH EXPAND comment before and after $mnh_expand$
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
# Count number of space blank before a first character of a line
    count=0
    for t in text:
        if t == " ":
            count=count+1
        else:
            break
    return count
        
def check_indent(n,indent_score,worktext):
# Correct the indentation with respect to the indent_score
    nblank=count_blank(worktext)
    rawtext=worktext[nblank:]
    correcttext=' '*((indent_score)*2) + rawtext
    return correcttext

def firstchar(string):
# Return the first non-blank character of a string
    for i in string:
        if i == " ": pass
        else: return i

def first7char(string):
# Return the 7 first non-blank character of a string
    nb_blank=0
    for i in string:
        if i == " ": nb_blank=nb_blank+1 
        else: return string[nb_blank:nb_blank+7]


def correct_indent(f):
    import sys
# Correct the indentation between MNH EXPAND comment
# Does not change the indentation outside MNH EXPAND > ... < MNH END EXPAND
# Does not indent comment lines
# Handles IF in two lines (e.g IF(... &
#                                 ......) THEN                                    
# TODO : handles more than two lines of & with IF
# TODO : do not correct indentation for & within an regular line
# To improve : the indentation correction is weird (strict) if the indentation before
# the MNH EXPAND comments is not respected (over-indentation)
    fin = open(f,'r')
    fout  = open(f+'_CORRECT_INDENT','w')
    contentbyline = fin.readlines()
    
    ncurrline=0
    indent_score=0
    expand_score=0
    passNextLine = {'Pass':False, 'Reason':"" }
    for i in contentbyline:
        i7=first7char(i)
        if passNextLine['Pass']: #Second line with & for if or #if(n)def
            textwrite=i
            passNextLine['Pass']=False
            if passNextLine['Reason'] == "IF":
                indent_score = indent_score+1
            elif passNextLine['Reason'] == "#ifdef":
                textwrite=check_indent(ncurrline,indent_score,i)
            elif passNextLine['Reason'] == "#else":
                textwrite=check_indent(ncurrline,indent_score,i)
                # ONLY IF present in between #ifdef is handled
                # If more test is needed (present in the fortran code), duplicate test here
                if "IF" in i7 and "THEN" in i:
                    indent_score = indent_score+1
            else:
                 sys.exit("Reason for passing the line not defined")
        elif "! $MNH EXPAND$ !" in i:
            expand_score = expand_score + 1
            textwrite = ""
        elif "! $MNH END EXPAND$ !" in i:
            expand_score = expand_score - 1
            textwrite = ""
        # Ignore comment lines
        elif firstchar(i) == "!":
            textwrite=i
        # Correct indentation only in between $MNH EXPAND$ and $MNH END EXPAND$
        elif expand_score >= 1:
            if indent_score > 0 and ("DO J" not in i and "ENDDO" not in i and "END DO" not in i and "END IF" not in i and "ENDIF" not in i and "THEN" not in i and "ELSE" not in i and "#if" not in i and "#else" not in i and "#endif" not in i):
                textwrite=check_indent(ncurrline,indent_score,i)
            # #ifdef handling = pass to next lines with no indent
            elif "#if" in i7:  #ifdef or ifndef
                print(i)
                textwrite=i
                passNextLine['Pass'],passNextLine['Reason']= (True, "#ifdef")
            elif "#else" in i7:
                textwrite=i
                passNextLine['Pass'],passNextLine['Reason']= (True, "#else")
            elif '#endif' in i7:
                textwrite=i
            elif "ELSE" in i7: #ELSE or ELSEIF
                indent_score = indent_score-1        
                textwrite=check_indent(ncurrline,indent_score,i)
                indent_score = indent_score+1
            elif "IF" in i7 and "&" in i: #IF on two lines #TODO on > 2 lines
                textwrite=i
                passNextLine['Pass'],passNextLine['Reason']= (True, "IF")
            elif "IF" in i7 and "THEN" in i: #exclude IF in one line (without THEN)
                textwrite=check_indent(ncurrline,indent_score,i)
                indent_score = indent_score+1        
            elif "DO J" in i7:
                textwrite=check_indent(ncurrline,indent_score,i)
                indent_score = indent_score+1
            elif ("ENDDO" in i7) or ("END DO" in i7):
                indent_score = indent_score-1
                textwrite=check_indent(ncurrline,indent_score,i)
            elif ("END IF" in i7) or ("ENDIF" in i7):
                indent_score = indent_score-1
                textwrite=check_indent(ncurrline,indent_score,i)       
            else:
                textwrite=i
        else: # not EXPAND lines nor within mnh_expand
            if "ELSE" in i7: #ELSE or ELSEIF
                pass #no increase in indent
            elif "IF" in i7 and "&" in i: #IF on two lines #TODO on > 2 lines
                passNextLine['Pass'],passNextLine['Reason']= (True, "IF")
            elif "IF" in i7 and "THEN" in i: #exclude IF in one line (without THEN)
                indent_score = indent_score+1        
            elif "DO J" in i7:
                indent_score = indent_score+1
            elif ("ENDDO" in i7) or ("END DO" in i7):
                indent_score = indent_score-1
            elif ("END IF" in i7) or ("ENDIF" in i7):
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
