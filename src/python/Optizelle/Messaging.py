# Defines how we output messages to the user"

import sys
import Optizelle.Exception

#---Messaging0---
def t(msg):
    """At it's core, we take in a string and then write it somewhere"""
    raise Optizelle.Exception.t("Messaging function not defined")
#---Messaging1---

def stdout(msg):
    """Write a string to stdout"""
    sys.stdout.write("%s\n" % msg)
