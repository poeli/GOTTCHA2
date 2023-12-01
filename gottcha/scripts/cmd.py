
#!/usr/bin/env python3
import sys
import io
from gottcha.scripts import pull_database
from gottcha.scripts import gottcha2

def usage():
    print("usage: gottcha2 profile --help\ngottcha2 pull --help")
    sys.exit()

def gottcha2_command():
    args = sys.argv[1:]
    if len(args) < 1:
        usage()
    elif args[0] == "profile":
        gottcha2.main(args[1:])
    elif args[0] == "pull":
        pull_database.main(args[1:])
    else:
        print("%s is not a valid argument", args[1])
        usage()
