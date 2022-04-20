package main

import (
	"fmt"
	"math"
	"strconv"
)

func problem1() {
	upperLimit := 1000
	var sum int
	for i := 0; i < upperLimit; i++ {
		if (i%3 == 0) || (i%5 == 0) {
			sum += i
		}
	}
	fmt.Println(sum)
}

func fibonacci(n int) int {
	if n <= 0 {
		return 1
	} else if n == 1 {
		return 2
	} else {
		return fibonacci(n-1) + fibonacci(n-2)
	}
}

func isPrime(n int) bool {
	var isComp bool
	for i := 2; float64(i) <= math.Sqrt(float64(n)); i++ {
		isComp = isComp || (n%i == 0)
	}
	return !isComp
}

func smallestPrime(n int) int {
	var out int
	if isPrime(n) {
		out = n
	} else {
		for i := 2; float64(i) <= math.Sqrt(float64(n)); i++ {
			if (n%i == 0) && isPrime(i) {
				out = i
				break
			}
		}
	}
	return out
}

func primeFactorize(n int) []int {
	var output []int
	running := n
	for !isPrime(running) {
		output = append(output, smallestPrime(running))
		running = running / smallestPrime(running)
	}
	output = append(output, running)
	return output
}

func slicesEqual(s, t []int) bool {
	out := true
	if len(s) != len(t) {
		out = false
	} else {
		for i := 0; i < len(s); i++ {
			out = out && (s[i] == t[i])
		}
	}
	return out
}

func reverse(s []int) []int {
	out := make([]int, len(s))
	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		out[i], out[j] = s[j], s[i]
	}
	return out
}

func digits(n int) []int {
	var outs []int
	for i := 0; math.Pow10(i) <= float64(n); i++ {
		outs = append(outs, int(math.Floor(float64(n)/math.Pow10(i)))%10)
	}
	return reverse(outs)
}

func isPalindrome(n int) bool {
	return slicesEqual(digits(n), reverse(digits(n)))
}

func problem2() {
	upperLimit := 4000000
	var sum int
	for i := 0; fibonacci(i) < upperLimit; i++ {
		if i%3 == 1 {
			sum += fibonacci(i)
		}
	}
	fmt.Println(sum)
}

func problem3() {
	N := 600851475143
	facs := primeFactorize(N)
	fmt.Println(facs[len(facs)-1])
}

func problem4() {
	var found bool
	for i := 0; i < 1998 && !found; i++ {
		for j := 0; j <= i; j++ {
			if isPalindrome((999 - j) * (999 - (i - j))) {
				found = true
				fmt.Println(999-j, 999-(i-j), (999-j)*(999-(i-j)))
				break
			}
		}
		if found {
			break
		}
	}
}

func binGCD(n, m int) int {
	var a, b int
	if n <= m {
		a, b = m, n
	} else {
		a, b = n, m
	}
	for a%b != 0 {
		a, b = b, a%b
	}
	return b
}

func binLCM(n, m int) int {
	return (n * m) / binGCD(n, m)
}

func lcm(s []int) int {
	out := 1
	for i := 0; i < len(s); i++ {
		if s[i] < 0 {
			out = binLCM(out, -s[i])
		} else if s[i] > 0 {
			out = binLCM(out, s[i])
		}
	}
	return out
}

func problem5() {
	N := 20
	var firstN []int
	for i := 1; i <= N; i++ {
		firstN = append(firstN, i)
	}
	fmt.Println(lcm(firstN))
}

func problem6() {
	N := 100
	fmt.Println(((N)*(N+1)/2)*((N)*(N+1)/2) - ((N)*(N+1)*(1+(2*N)))/6)
}

func problem7() {
	N := 10001
	var primes []int
	for i := 2; len(primes) < N; i++ {
		if isPrime(i) {
			primes = append(primes, i)
		}
	}
	fmt.Println(primes[N-1])
}

func problem8() {
	bigNum := "7316717653133062491922511967442657474235534919493496983520312774506326239578318016984801869478851843858615607891129494954595017379583319528532088055111254069874715852386305071569329096329522744304355766896648950445244523161731856403098711121722383113622298934233803081353362766142828064444866452387493035890729629049156044077239071381051585930796086670172427121883998797908792274921901699720888093776657273330010533678812202354218097512545405947522435258490771167055601360483958644670632441572215539753697817977846174064955149290862569321978468622482839722413756570560574902614079729686524145351004748216637048440319989000889524345065854122758866688116427171479924442928230863465674813919123162824586178664583591245665294765456828489128831426076900422421902267105562632111110937054421750694165896040807198403850962455444362981230987879927244284909188845801561660979191338754992005240636899125607176060588611646710940507754100225698315520005593572972571636269561882670428252483600823257530420752963450"
	var numDigs []int
	for i := 0; i < len(bigNum); i++ {
		dig, _ := strconv.Atoi(bigNum[i : i+1])
		numDigs = append(numDigs, dig)
	}
	N := 13
	var prod, maxProd int
	for i := 0; i+N < len(bigNum); i++ {
		prod = 1
		seg := numDigs[i : i+N]
		for j := 0; j < N; j++ {
			prod *= seg[j]
		}
		if maxProd < prod {
			maxProd = prod
		}
	}

	fmt.Println(maxProd)
}

func main() {
	problem8()
}
