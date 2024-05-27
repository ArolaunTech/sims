stage0 = 1;
stage1 = 0;
stage2 = 0;
stage3 = 0;
stage4 = 0;
stage5 = 0;
stage6 = 0;
fullGrown = 0;
Ticks = 0;
fullGrowneff = 0;

while True:
	for i in range(20):
		fullGrown += 0.0205078125 * stage6
		stage6 += 0.0205078125 * (stage5 - stage6)
		stage5 += 0.0205078125 * (stage4 - stage5)
		stage4 += 0.0205078125 * (stage3 - stage4)
		stage3 += 0.0205078125 * (stage2 - stage3)
		stage2 += 0.0205078125 * (stage1 - stage2)
		stage1 += 0.0205078125 * (stage0 - stage1)
		stage0 *= 0.9794921875
		Ticks += 1
	print(Ticks)
	print(fullGrown)
	fullGrowneff = fullGrown / Ticks
	print(fullGrowneff)
	h = input("Continue?")
	if h.lower() == "yes":
		continue
	else:
		exit(0)

#What is supposed to be a simulator. The editor is awesome!
