def readHeader(inFile,EOH='EOH'):
    """Read and return the header of a text file

    If no header is present, the entire file is returned
    
    Input: the name of a newly opened text file
    EOH: acceptable 'end of header' lines
    
    Output:
    1) a list of strings, one per header line
    2) a boolean that is True if an EOH was found

    Intended uses:
    1) Read and print the header of a file.
    2) Read past the header in preparation for reading the data.
    3) Determine the absence of a header.

    Note: Expects newline character (octal 012 = \n).
    A supplemental carriage return is OK (octal 014 = \m or \r).
    Carriage return alone is not (Mac Excel does this when asked to produce
    CSV or MSDOSCSV, but not WindowsCSV.)
    
    written 6/2009 by LAM
    revised 2/26/2012 by LAM to add second EOH string and to output boolean
    revise 9/1/2019 by LAM to not use variable name 'file'
    revised 7/3/204 by LAM to allow an alternative EOH to be specified
    """
    header = list()
    while True:
        line = inFile.readline()
        if line[0:len(EOH)] == EOH: # EOH reached, header returned
            return(header, True)
        else:
            if line == "":                  # EOF reached, entire file returned
                return(header, False)
            else:                           # Another header line read
                header.append(line)

if __name__ == '__main__':
    #Run "readHeader", print results to standard output
    #User input: complete file name from which to read
    fileName = '../PythonData/Becker/BeckerPeriods.csv'
    #fileName = '../PythonData/Test/CalCSS-J0811-model.csv'

    #Open file; read header; print header
    inFile = open(fileName)
    (header, HasEOH) = readHeader(inFile)
    for line in header:
        print(line.strip('\n'))
    if HasEOH:
        print('EOH found')
    else:
        print('No EOH found')   
    print('\nreadHeader finished')
