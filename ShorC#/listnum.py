def factors(n):
    while n > 1:
        for i in range(2, n + 1):
            if n % i == 0:
                n //= i
                yield i
                break

s = ''

for n in range(64):
	l = list(factors(n))
	if not n%2 or len(l) != 2:
		continue
	s += str(n) + ','

print(s)