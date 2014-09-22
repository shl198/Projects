import subprocess
def Message(string,email):
    """
    This function send message to email when it run. 
    Used to calculate the time code runs.
    """
    cmd = ('echo {quote}|mailx -s "{string}" {email}').format(quote="",string=string,email=email)
    subprocess.call(cmd,shell=True)
