import random

def generate_coeffs_answer_n(n = 0):
	if n <=0:
		n = random.randrange(100+1)
	sol =  []
	for i in range(n):
		sol.append(random.randrange(-100,100))
#		if i == 0:
#			sol[0] = 0
		print(sol[-1])
	c = []
	a = []
	for i in range(n):
		coeffs = []
		for j in range(n):
			coeffs.append(random.randrange(-10,10))
#			if i == 0 and j ==  0:
#				coeffs[0] = 0
		answer = 0
		for j in range(n):
			answer += coeffs[j]*sol[j]
		c.append(coeffs)
		a.append(answer)
	return c,a,n

def write_trial(c,a,n, dest = None):
	lines = []
	lines.append(str(n))
	for i in range(n):
		temp = []
		for j in c[i]:
			temp.append(str(j))
		lines.append(' '.join(temp))
	for i in range(n):
		lines.append(str(a[i]))
	if dest == None:
		for i in lines:
			print(i)
	else:
		with open(dest, 'w') as f:
			for i in lines:
				f.write(str(i) + "\n")
def main():
	c,a,n = generate_coeffs_answer_n(4)
	write_trial(c,a,n,"data.dat")
main()
