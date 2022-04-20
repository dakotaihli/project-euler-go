package main

import (
	"fmt"
	"math"
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

func main() {
	problem6()
}
