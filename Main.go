package main

import (
	"errors"
	"math/big"
	"sort"
	"strconv"
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

func intPow(x, n int) int {
	if n <= 0 {
		return 1
	} else if n == 1 {
		return x
	} else {
		return x * intPow(x, n-1)
	}
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

// isPrime returns whether |n| is prime
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

// primeFactorize returns the prime factors of n
// If withDuplicates is false, only return distinctfactors
func primeFactorize(n int, withDuplicates bool) []int {
	var output []int
	running := n
	for running != 1 {
		if withDuplicates || len(output) == 0 || output[len(output)-1] != smallestPrime(running) {
			output = append(output, smallestPrime(running)) // If withDuplicates is false, this will only add this prime once
		}
		running /= smallestPrime(running)
	}
	return output
}

func primeAfter(n int) int {
	for p := n + 1; p <= 2*n; p++ {
		if isPrime(p) {
			return p
		}
	}
	return 0 //This will never happen since there is always a prime between n and 2n
}

func distinctPrimeFactors(n int) []int {
	var output []int
	running := n
	for running != 1 {
		p := smallestPrime(running)
		output = append(output, p)
		for running%p == 0 {
			running /= p
		}
	}
	return output
}

// primesBelow returns a sorted slice of ints, consisting of all primes strictly below the input
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
	var out bool
	if len(s) != len(t) {
		out = false
	} else {
		out = true
		for i := 0; i < len(s); i++ {
			out = out && (s[i] == t[i])
		}
	}
	return out
}

func sortedIntSliceContains(s []int, x int) bool {
	return x == s[sort.SearchInts(s, x)]
}

func reverse(s string) string {
	t := []rune(s)
	for i, j := 0, len(t)-1; i <= j; i, j = i+1, j-1 {
		t[i], t[j] = t[j], t[i]
	}
	return string(t)
}

func reverseInt(s []int) []int {
	out := make([]int, len(s))
	for i, j := 0, len(s)-1; i <= j; i, j = i+1, j-1 {
		out[i], out[j] = s[j], s[i]
	}
	return out
}

func reverseString(s []string) []string {
	out := make([]string, len(s))
	for i, j := 0, len(s)-1; i <= j; i, j = i+1, j-1 {
		out[i], out[j] = s[j], s[i]
	}
	return out
}

func numToDigits(n int) []int {
	var outs []int
	nRunning := n
	for i := 0; intPow(10, i) <= n; i++ {
		outs = append(outs, nRunning%10)
		nRunning /= 10
	}
	return reverseInt(outs)
}

func digitsToNum(s []int) int {
	power := 1
	var num int
	for i := 0; i < len(s); i++ {
		num += s[len(s)-1-i] * power
		power *= 10
	}
	return num
}

func areNumsPerms(n, m int) bool {
	nDigs, mDigs := numToDigits(n), numToDigits(m)
	sort.Ints(nDigs)
	sort.Ints(mDigs)
	return slicesEqual(nDigs, mDigs)
}

func removeSliceElt(s []int, i int) []int {
	if (i < 0) || i >= len(s) {
		return s
	} else {
		out := make([]int, len(s)-1)
		for j := 0; j < len(out); j++ {
			if j < i {
				out[j] = s[j]
			} else {
				out[j] = s[j+1]
			}
		}
		return out
	}
}

func intPerms(s []int) [][]int {
	var out [][]int
	if len(s) == 0 {
		return [][]int{nil}
	} else {
		for i, x := range s {
			sLess := removeSliceElt(s, i)
			t := intPerms(sLess)
			for _, y := range t {
				out = append(out, append(y, x))
			}
		}
		return out
	}
}

func intSliceDistinctElts(s []int) []int {
	var out []int
	for _, x := range s {
		var contains bool
		for _, y := range out {
			contains = contains || (x == y)
		}
		if !contains {
			out = append(out, x)
		}
	}
	return out
}

func isPalindromeInt(n int) bool {
	return slicesEqual(numToDigits(n), reverseInt(numToDigits(n)))
}

func isPalindromeString(s string) bool {
	isPal := true
	for i, j := 0, len(s)-1; i < j; i, j = i+1, j-1 {
		isPal = isPal && s[i] == s[j]
	}
	return isPal
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

func nCr(n, r int) int {
	if r < 0 || r > n || n < 0 {
		return 0
	} else if r == 0 || r == n {
		return 1
	} else {
		return nCr(n-1, r) + nCr(n-1, r-1)
	}
}

func triangle(n int) int {
	return (n * (n + 1)) / 2
}

func pentagonal(n int) int {
	return (n * (3*n - 1)) / 2
}

func hexagonal(n int) int {
	return triangle(2*n - 1)
}

func pairIterate(x, y int) (int, int) {
	if x == 0 {
		return y + 1, 0
	} else {
		return x - 1, y + 1
	}
}

// Used for Problems 18 and 67
func trianglePathReduce(tri [][]int) ([][]int, error) {
	for i, l := range tri {
		if len(l) != i+1 {
			return tri, errors.New("triangle is incorrectly formatted")
		}
	}
	var newTri [][]int = tri[:len(tri)-2]
	newLev := make([]int, len(tri)-1)
	for i := 0; i < len(newLev); i++ {
		newLev[i] = tri[len(tri)-2][i] + max(tri[len(tri)-1][i], tri[len(tri)-1][i+1])
	}
	return append(newTri, newLev), nil
}

/*
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
*/

func collatz(n int) int {
	if n%2 == 0 {
		return n / 2
	} else {
		return (3 * n) + 1
	}
}

func digitSum(s string) int {
	var sum int
	for i := 0; i < len(s); i++ {
		d, _ := strconv.Atoi(s[i : i+1])
		sum += d
	}
	return sum
}

// Only works for 1 <= n <= 1000
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

func bigFactOnBig(z *big.Int) *big.Int {
	if z.Cmp(big.NewInt(int64(0))) == -1 {
		return big.NewInt(int64(0))
	} else {
		running := big.NewInt(int64(1))
		for i := big.NewInt(int64(1)); i.Cmp(z) < 1; i.Add(i, big.NewInt(int64(1))) {
			running.Mul(running, i)
		}
		return running
	}
}

func bigFact(n int) *big.Int {
	return bigFactOnBig(big.NewInt(int64(n)))
}

func stringToASCII(s string) []int {
	out := make([]int, len(s))
	for i, c := range s {
		out[i] = int(c)
	}
	return out
}

// alphValue is used for problems 22 and 42
func alphValue(s string) int {
	var sum int
	for _, l := range s {
		sum += (int(l) - 64)
	}
	return sum
}

func sumOfDivisors(n int) int {
	sum := 1
	for i := 2; i < (n/2)+1; i++ {
		if n%i == 0 {
			sum += i
		}
	}
	return sum
}

func isAbundant(n int) bool {
	return sumOfDivisors(n) > n
}

func isLeap(y int) bool {
	if y%4 != 0 {
		return false
	} else if y%100 == 0 && y%400 != 0 {
		return false
	} else {
		return true
	}
}

func isDateBeforeOrSame(ya, ma, da, yb, mb, db int) bool {
	return (ya < yb) || (ya == yb && ma < mb) || (ya == yb && ma == mb && da <= db)
}

func tomorrow(y, m, d int) (int, int, int) {
	if m == 2 {
		if (isLeap(y) && d == 29) || (!isLeap(y) && d == 28) {
			return y, 3, 1
		} else {
			return y, m, d + 1
		}
	} else if m == 4 || m == 6 || m == 9 || m == 11 {
		if d == 30 {
			return y, m + 1, 1
		} else {
			return y, m, d + 1
		}
	} else if m == 12 && d == 31 {
		return y + 1, 1, 1
	} else {
		if d == 31 {
			return y, m + 1, 1
		} else {
			return y, m, d + 1
		}
	}
}

func allPos(in []int) bool {
	out := true
	for _, n := range in {
		out = out && (n > 0)
	}
	return out
}

func coinCombos(coins []int, total int) [][]int {
	if allPos(coins) && total >= 0 {
		if total == 0 {
			emptySoln := make([]int, len(coins))
			return [][]int{emptySoln}
		} else if len(coins) == 0 {
			return make([][]int, 0)
		} else {
			outs := make([][]int, 0)
			for i := 0; i*coins[len(coins)-1] <= total; i++ {
				partialSols := coinCombos(coins[:len(coins)-1], total-(i*coins[len(coins)-1]))
				for _, s := range partialSols {
					outs = append(outs, append(s, i))
				}
			}
			return outs
		}
	} else {
		return nil
	}
}

func isAdmissible(n int) bool {
	// One checks that n is admissible iff it is at least 2, and its distinct
	// prime factors form an initial segment of the primes
	if n%2 != 0 || n < 2 {
		return false
	} else {
		running := n
		for p := 2; running%p == 0; p = primeAfter(p) {
			for running%p == 0 {
				running /= p
			}
		}
		return running == 1
	}
}

func main() {
	problem(58)
}
