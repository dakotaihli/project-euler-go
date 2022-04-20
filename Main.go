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

func main() {
	problem3()
}
