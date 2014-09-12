import subprocess
def Message(string):
    """
    This function send message to email when it run. Used to calculate the time code runs.
    """
    cmd = ('echo {quote}|mailx -s "{string}" shl198@eng.ucsd.edu').format(quote="",string=string)
    subprocess.call(cmd,shell=True)
