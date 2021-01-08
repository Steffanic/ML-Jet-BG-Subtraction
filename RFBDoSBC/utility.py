def msg(_msg):
    print("JEB: %s"%_msg)

def log(_msg, filename="info.txt"):
    with open(filename, "a") as f:
        f.write(_msg+'\n')
