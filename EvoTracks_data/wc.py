def readData(file,delimiter=None):
    """Read and return the data from a file, then close the file.
    
    Input:
    1) the file name
    2) delimiter for parsing the line.  The default delimiter is white space.
    
    Output:
    1) a list of lists, each of which has the data from one line.

    Intended uses:
    1) Read in the data from a text file.  Typically, this will be csv data
    and the header will already have been passed over.  However, any delimiter
    may be used, and the status of the header is up to the user.
    2) wc.py uses this to report on size statistics.
    Note: a blank final line is not counted in the tally (so that the lack of
    a trailing newline character is  unimportant).
    
    written 7/2009 by LAM
    revised 7/6/2009 to use readHeader function 
    revised 4/21/2010 to strip the \n off the last item
    revised 2/26/2012 to use revised readHeader (which reports on EOH) and to
        give statistics on full file in absence of a header
    """
    data = list()
    line = file.readline()
    while line:
        if(line != ''):
            if delimiter == None:
                data.append(line.strip('\n').split())
            else:
                data.append(line.strip('\n').split(delimiter))
            line = file.readline()
    file.close()
    return(data)

if __name__ == '__main__':
    from readHeader import readHeader
    #Run "readData", print statistics on file contents akin to Linux wc command
    #User input: complete file name from which to read; delimiter
    fileName = '../PythonData/Test/CalCSS-J0811-model.csv'
    delimiter = ',' # appropriate for a comma separated file

    #Open file; read over header; read out data; report statistics
    file = open(fileName)
    (header, HasEOH) = readHeader(file)
    if not HasEOH:                  # if no EOH, close file and read from top
        file.close()
        file = open(fileName)
    data = readData(file, delimiter)
    if not HasEOH:
        print('No header found, statistics of entire file:')
    else:
        print('Header found, statistics of data only:')
    print('Lines of data:', len(data))
    count = 0
    for line in data:
        count += len(line)
    print('Data cells:', count)
    print('\nwc finished')
