import ftplib

def ftp_check(ftp, host, user, password, pwd):
    '''
    Pings the FTP server to make sure the connection is live,
    reconnects if it isn't.
    '''
    try:
        #ping server
        ftp.voidcmd("NOOP")
        return(ftp)
    #if connection has timed out
    except ftplib.error_temp:
        #reconnect
        ftp = ftp_connect(host, user, password, directory = pwd)
        return(ftp)

def ftp_connect(host, user, password, directory = None):
    '''
    Connect to FTP server.
    directory: if specified, change to that directory.
    '''
    connected = False
    while not connected:
        try:
            ftp = ftplib.FTP(host, timeout = 10000)
            connected = True
        except TimeoutError:
            print("TimeoutError! Trying again...")
    ftp.login(user, password)
    if directory:
        ftp.cwd(directory)
    return(ftp)

def ftp_retrieve(ftp, host, user, password, directory, file_name, destination = None):
    '''
    Retrieve one or several files from an FTP site.
    Meant to be given a live FTP connection, with the correct working directory, but still needs information to connect in case there is a timeout.
    directory: source directory on the FTP site (only used in case of timeout)
    file: name of file to retrieve
    destination: save the file to this location. If unspecified, the current working directory will be used.
    '''
    if destination:
        #this is to make it easier to join the directory path with a file name
        destination = "{0}/".format(destination)
    else:
        destination = ""
    local_file_name = "{0}{1}".format(destination, file_name)
    #it's this complicated because you want to be able to retrieve binary data
    with open(local_file_name, "wb") as local_file:
        #check that the connection is live, reconnect otherwise
        ftp = ftp_check(ftp, host, user, password, directory)
        retrieved = False
        #sometimes the file doesn't transfer properly so you have to keep on
        #trying till you get it
        while not retrieved:
            try:
                ftp.retrbinary("RETR {0}".format(file_name), local_file.write)
                retrieved = True
            except EOFError:
                print("EOFError! Trying again...")
                pass
            except TimeoutError:
                print("TimeoutError! Trying again...")
                ftp = ftp_check(ftp, host, user, password, directory)
    print("Retrieved file {0}.".format(file_name))
    return(ftp)
