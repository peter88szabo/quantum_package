#!/usr/bin/env python2

import sys

with open(sys.argv[1],'r') as f:
  raw_data = f.read()

output = []
inside = False
level = 0
for i in raw_data:
  new_i = i
  if i == "@":
    inside = True
  elif i == "{" and inside and level == 0:
    new_i = ""
  elif i == "}" and inside and level == 1:
    inside = False
    new_i = ""
  if i == "{":
    level += 1
  elif i == "}":
    level -= 1
  output.append(new_i)

print "".join(output).replace("@test","echo").replace("|| skip","|| return")



