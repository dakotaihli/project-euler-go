package main

import (
	"math"
	"strconv"
	"strings"
)

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func max(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func maxS(nums []int) int {
	maxVal := nums[0]
	for _, n := range nums {
		if n > maxVal {
			maxVal = n
		}
	}
	return maxVal
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

//isPrime returns whether |n| is prime
func isPrime(n int) bool {
	if n < 0 {
		n *= -1
	}
	if (n < 2) || (n > 3 && n%6 != 1 && n%6 != 5) {
		return false
	} else if n == 2 || n == 3 {
		return true
	} else {
		var isComp bool
		for i := 2; i*i <= n; i++ {
			isComp = isComp || (n%i == 0)
		}
		return !isComp
	}
}

func smallestPrime(n int) int {
	var out int
	if isPrime(n) {
		out = n
	} else {
		for i := 2; i*i <= n; i++ {
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

func primesBelow(n int) []int {
	var output []int
	for i := 2; i < n; i++ {
		if isPrime(i) {
			output = append(output, i)
		}
	}
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

func reverseInt(s []int) []int {
	out := make([]int, len(s))
	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		out[i], out[j] = s[j], s[i]
	}
	return out
}

func reverseString(s []string) []string {
	out := make([]string, len(s))
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
	return reverseInt(outs)
}

func isPalindrome(n int) bool {
	return slicesEqual(digits(n), reverseInt(digits(n)))
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

func divisors(n int) []int {
	var divs []int
	for i := 1; i <= n; i++ {
		if n%i == 0 {
			divs = append(divs, i)
		}
	}
	return divs
}

func triangle(n int) int {
	return (n * (n + 1)) / 2
}

//Need to add error handling if inputs aren't numbers
func stringAdd(summands []string) string {
	var reverseSum []int
	summandSlices := make([][]string, len(summands))
	lengths := make([]int, len(summands))
	for i, s := range summands {
		lengths[i] = len(s)
		summandSlices[i] = strings.Split(s, "")
	}
	maxLength := maxS(lengths)
	var carry int
	for i := 0; i < maxLength || carry > 0; i++ {
		sum := 0
		sum += carry
		for _, s := range summandSlices {
			if i < len(s) {
				d, _ := strconv.Atoi(s[len(s)-(i+1)])
				sum += d
			}
		}
		reverseSum = append(reverseSum, sum%10)
		carry = sum / 10
	}
	var out string
	for i := 0; i < len(reverseSum); i++ {
		out = out + strconv.Itoa(reverseSum[len(reverseSum)-(i+1)])
	}
	return out
}

func binStringMult(s, t string) string {
	if len(s) < len(t) {
		s, t = t, s
	}
	tSlice := strings.Split(t, "")
	var auxProds []string
	for i := 0; i < len(t); i++ {
		jDig, _ := strconv.Atoi(tSlice[len(t)-(i+1)])
		sTimesPow := s
		for z := 0; z < i; z++ {
			sTimesPow = sTimesPow + "0"
		}
		for j := 0; j < jDig; j++ {
			auxProds = append(auxProds, sTimesPow)
		}
	}
	return stringAdd(auxProds)
}

func stringMult(nums []string) string {
	if len(nums) == 0 {
		return "1"
	} else {
		running := nums[0]
		for i := 1; i < len(nums); i++ {
			running = binStringMult(running, nums[i])
		}
		return running
	}
}

func stringPow(base string, pow int) string {
	var bases []string
	for i := 1; i <= pow; i++ {
		bases = append(bases, base)
	}
	return stringMult(bases)
}

func collatz(n int) int {
	if n%2 == 0 {
		return n / 2
	} else {
		return (3 * n) + 1
	}
}

func digitSum(s string) string {
	return stringAdd(strings.Split(s, ""))
}

//Only works for 1 <= n <= 1000
func numberName(n int) string {
	digits := []string{"one", "two", "three", "four", "five", "six", "seven", "eight", "nine"}
	teens := []string{"ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen", "sixteen", "seventeen", "eighteen", "nineteen"}
	tennies := []string{"twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty", "ninety"}
	if n == 1000 {
		return "one thousand"
	} else {
		var out string
		nModHund := n % 100
		if n >= 100 {
			out = digits[(n/100)-1] + " hundred"
			if nModHund != 0 {
				out = out + " and "
			}
		}
		if nModHund/10 == 0 && nModHund != 0 {
			out = out + digits[(n%10)-1]
		} else if nModHund/10 == 1 {
			out = out + teens[n%10]
		} else if nModHund/10 > 1 {
			out = out + tennies[(nModHund/10)-2]
			if n%10 != 0 {
				out = out + "-" + digits[(n%10)-1]
			}
		}
		return out
	}
}

func main() {
	problem(48)
}
