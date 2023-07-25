def file_cmp():
    # reading files
    f1 = open("/home/anurag/PycharmProjects/Ibab/output_files/mutated_ligplot.sum", "r")
    f2 = open("/home/anurag/PycharmProjects/Ibab/output_files/altered_ligplot.sum", "r")
    f3 = open("/home/anurag/PycharmProjects/Ibab/ligplot_differences.txt", "w")

    f1_data = f1.readlines()
    f2_data = f2.readlines()
    syntax = f1_data[1:5]

    f3.write("Interactions Gained: \n")
    f3.writelines(syntax)

    key=0

    for line1 in f1_data[1:]:
        if line1 not in f2_data:
            key=1
            f3.write(line1)
    if key == 0:
        f3.write("None \n")

    f3.write("\n")
    f3.write("\n")

    key=0
    f3.write("Interactions Lost: \n ")
    f3.writelines(syntax)

    for line2 in f2_data[1:]:
        if line2 not in f1_data:
            key=1
            f3.write(line2)
    if key == 0:
        f3.write("None \n")

    # closing files
    f1.close()
    f2.close()
    f3.close()

file_cmp()
