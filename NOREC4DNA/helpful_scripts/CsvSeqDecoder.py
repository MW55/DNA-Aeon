import pandas

csv = pandas.read_csv("/home/michael/Code/norec4dna/unilogo/results/finalData/filtered_table.csv")

is_bigger_th = csv #[csv["UniLogo_A"] < 100]
print(is_bigger_th)
res = is_bigger_th[is_bigger_th["sequences"].str.len() == 164]
print(res["sequences"])
i = 0
for x in res["sequences"]:
    i += 1
    with open("out/" + str(i) + ".RU10_DNA", "w") as f:
        f.write(x.upper())