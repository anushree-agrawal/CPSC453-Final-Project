file = "lunggenes"
text = open(file, 'r')
writefile = open("output.txt", 'w')
writefile.write("TEXT.PnS={")

for line in text:
    if line[-1] == '\n':
        line = line[:-1]
    writefile.write("'")
    writefile.write(line)
    writefile.write("',")
writefile.write("}")
