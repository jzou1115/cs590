labels= open("labels2.txt", "r").readlines()
results= open("results.txt", "r").readlines()

results2=open("results2.txt", "w")
for r in results:
    results2.write(labels[int(r)])
results2.close()



