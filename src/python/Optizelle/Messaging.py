# Defines how we output messages to the user"

import sys
import Optizelle.Exception

#---Messaging0---
def t(msg):
    """At its core, we take in a string and then write it somewhere"""
    raise Optizelle.Exception.t("Undefined messaging function")
#---Messaging1---

def stdout(msg):
    """Write a string to stdout"""
    print msg
