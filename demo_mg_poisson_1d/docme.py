import pydoc,re, sys, os, glob, signal

def bold(text):
    """Format a string in bold."""
    return "\033[1m" + text +"\x1b[0;0m"

def italic(text):
    """Format a string in italic."""
    return '\x1b[3m' + text + '\x1b[0m'

def underline(text):
    """Format a string in underline."""
    return "\x1b[4m"+text+"\x1b[0m"

def docme(path):
    '''Eextract the information for documentation of the file.'''
    text = ''
    for filename in glob.glob(os.path.join(path, '*.c')):
        src = open(filename, "r") #open source file
        text += ''.join(src.readlines())
        src.close()
        #pydoc.pager(text)
    
    #italic and color the text
    #print('\x1b[3;31;43m' + 'Hello world!' + '\x1b[0m') 
    comment = re.compile(r'\/\*(?P<comment>(?:[^*]+|\*[^/])+)\*\/')
    
    #getparbool/getparlargeint/getparint/getparfloat/getpardouble
    param = re.compile(r'(?:if\s*\(\!)?\s*getpar'
                       '(?P<type>largeint|int|float|double)'
                       '\s*\(\s*\"(?P<name>\w+)\"\s*\,'
                       '\s*\&(?P<var>[\w\_\[\]\.\+\-\>]+)\s*[\)]\s*[\)]?\s*'
                       '(?:[\{]|' # either \{ or
                       '(?:(?P=var)\s*\=\s*(?P<default>[^\;]+)|'
                       'err[^\;]+)?' # or err
                       '[\;])\s*' # ending with ;
                       '(?:\/\*\s*(?P<range>[\[][^\]]+[\]])?\s*'
                       '(?P<desc>(?:[^*]|\*[^/])+)\*\/)?') # comment
    #getparstring
    stringpar =  re.compile(r'(?:if\s*\(\!)?\s*getparstring'
                            '\s*\(\s*\"(?P<name>\w+)\"\s*\,'
                            '\s*\&(?P<var>\w+)\s*'
                            '(?:[\{]|' # either \{ or
                            '(?:(?P=var)\s*\=\s*(?P<default>[^\;]+)|'
                            'err[^\;]+)?' # or err
                            '[\;])\s*' # ending with ;
                            '(?:\/\*\s*(?P<range>[\[][^\]]+[\]])?\s*'
                            '(?P<desc>(?:[^*]|\*[^/])+)\*\/)?') # comment

    synopsis = re.compile(r'\s*Takes\s*\:\s*((?:[^\n]|[\n][^\n])+)'
                          '((?:.|\n)*)$')
    
    #---------------------------------------------------------------
    print bold('name').upper()
    name = '\t'
    print name
    
    #---------------------------------------------------------------
    print bold('description').upper()
    description = '\t'
    print description
    
        
    #---------------------------------------------------------------
    print bold("synopsis").upper()
    parline = '\t'
    pars = param.findall(text)
    for par in pars:
        type = par[0]
        parname = par[1]
        default = par[3]
        range = par[4]
        desc = par[5]
        parline += " %s=%s" % (parname,default)
    
    pars = stringpar.findall(text)
    for par in pars:
        type = 'string'
        parname = par[0]
        desc = par[1]
        parline += " %s=%s" % (parname,default)
        
    print parline

    #---------------------------------------------------------------
    print bold('parameters').upper()
    
    pars = param.findall(text)
    for par in pars:
        type = par[0]
        parname = par[1]
        default = par[3]
        range = par[4]
        desc = par[5]
        
        print '\t',underline(type),'\t', bold(parname+'='+default),range,desc
        
    pars = stringpar.findall(text)
    for par in pars:
        type = 'string'
        parname = par[1]
        default = par[3]
        range = par[4]
        desc = par[5]
        
        print '\t',underline(type),'\t',bold(parname+'='),desc


    #---------------------------------------------------------------
    print bold("see also").upper()
    seealso = '\t'
    print seealso

    #---------------------------------------------------------------
    print bold("author").upper()
    print '\t'+'Pengliang Yang, 2019, Email: ypl.2100@gmail.com '

    
if __name__ == '__main__':
    path = '.'
    docme(path)


