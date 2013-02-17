#!/usr/bin/python
#montage -background black +frame +shadow +label -tile 5x1 -geometry 200x300+0+0 A-done.png B-done.png C-done.png D-done.png E-done.png joined.png

import subprocess
import cgi
import base64
import string

form = cgi.FieldStorage()
if "value" not in form:
	print "Location: /valentines-2013/"
	print
	exit()

print "Content-Type: text/html"
print
#print form.getvalue("value")

cmd = "/usr/bin/montage -background black +frame +shadow +label -tile %ix1 -geometry 200x300+0+0"%len(form.getvalue("value"))
short_str = ""
for c in form.getvalue("value"):
	if c == " ":
		cmd += " /home/apache/p/htdocs/valentines-2013/space.png"
	elif c in string.ascii_letters:
		cmd += " /home/apache/p/htdocs/valentines-2013/%s-done.png"%c.upper()
		short_str += c.upper()
import random
v = random.randint(1,1000000)
cmd += " /home/apache/p/htdocs/valentines-2013/tmp/valentines-%s-%s.png"%(short_str,v)

subprocess.call(cmd, shell=True)

#print open("/home/apache/p/htdocs/valentines-2013/tmp/valentines-%s.png"%v, "rb").read()
print """
<html><head>
<link rel="stylesheet" href="/cloud/static/css/screen.css" type="text/css" media="screen, projection">
<link rel="stylesheet" href="/cloud/static/css/print.css" type="text/css" media="print">
</head>"""
print "<body style='background:black'><div style=''><img style='margin-left:auto;margin-right:auto' align='middle' src='/valentines-2013/tmp/valentines-%s-%s.png'/><br/><a href='/valentines-2013/'>Make your own...</a></div>"%(short_str,v)
print """
<div id="footer" style="text-align:left;position:absolute;bottom:0;">
<small style="margin-bottom:5px"><i><a href="http://pillow.rscheme.org">Proudly Pillow Powered</a>. Copyright(c) 2013 by Lane Kolbly</i></small>
</div>
"""
print "</body></html>"
