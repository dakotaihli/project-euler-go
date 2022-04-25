package main

import (
	"fmt"
	"math"
	"strconv"
	"strings"
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

func problem9() {
	var ag, bg, cg int
	for c := 334; c < 998; c++ {
		for b := 2; b+c < 999; b++ {
			if (b*b)+(1000-(b+c))*(1000-(b+c)) == c*c {
				ag, bg, cg = (1000 - (b + c)), b, c
			}
		}
	}
	fmt.Println(ag * bg * cg)
}

func problem10() {
	n := 2000000
	primes := primesBelow(n)
	var sum int
	for _, p := range primes {
		sum += p
	}
	fmt.Println(sum)
}

func problem11() {
	tableString := "08 02 22 97 38 15 00 40 00 75 04 05 07 78 52 12 50 77 91 08\n49 49 99 40 17 81 18 57 60 87 17 40 98 43 69 48 04 56 62 00\n81 49 31 73 55 79 14 29 93 71 40 67 53 88 30 03 49 13 36 65\n52 70 95 23 04 60 11 42 69 24 68 56 01 32 56 71 37 02 36 91\n22 31 16 71 51 67 63 89 41 92 36 54 22 40 40 28 66 33 13 80\n24 47 32 60 99 03 45 02 44 75 33 53 78 36 84 20 35 17 12 50\n32 98 81 28 64 23 67 10 26 38 40 67 59 54 70 66 18 38 64 70\n67 26 20 68 02 62 12 20 95 63 94 39 63 08 40 91 66 49 94 21\n24 55 58 05 66 73 99 26 97 17 78 78 96 83 14 88 34 89 63 72\n21 36 23 09 75 00 76 44 20 45 35 14 00 61 33 97 34 31 33 95\n78 17 53 28 22 75 31 67 15 94 03 80 04 62 16 14 09 53 56 92\n16 39 05 42 96 35 31 47 55 58 88 24 00 17 54 24 36 29 85 57\n86 56 00 48 35 71 89 07 05 44 44 37 44 60 21 58 51 54 17 58\n19 80 81 68 05 94 47 69 28 73 92 13 86 52 17 77 04 89 55 40\n04 52 08 83 97 35 99 16 07 97 57 32 16 26 26 79 33 27 98 66\n88 36 68 87 57 62 20 72 03 46 33 67 46 55 12 32 63 93 53 69\n04 42 16 73 38 25 39 11 24 94 72 18 08 46 29 32 40 62 76 36\n20 69 36 41 72 30 23 88 34 62 99 69 82 67 59 85 74 04 36 16\n20 73 35 29 78 31 90 01 74 31 49 71 48 86 81 16 23 57 05 54\n01 70 54 71 83 51 54 69 16 92 33 48 61 43 52 01 89 19 67 48"
	tableRows := strings.Split(tableString, "\n")
	var tableCells [][]string
	for i := 0; i < len(tableRows); i++ {
		tableCells = append(tableCells, strings.Split(tableRows[i], " "))
	}
	tableNums := make([][]int, len(tableCells))
	for i := 0; i < len(tableCells); i++ {
		tableNums[i] = make([]int, len(tableCells[i]))
		for j := 0; j < len(tableCells[i]); j++ {
			tableNums[i][j], _ = strconv.Atoi(tableCells[i][j])
		}
	}
	var biggestHoriz, biggestVert, biggestDiag, biggestRDiag int
	for i := 0; i < len(tableNums); i++ {
		for j := 0; j < len(tableNums[i]); j++ {
			if i+3 < len(tableNums) && tableNums[i][j]*tableNums[i+1][j]*tableNums[i+2][j]*tableNums[i+3][j] > biggestVert {
				biggestVert = tableNums[i][j] * tableNums[i+1][j] * tableNums[i+2][j] * tableNums[i+3][j]
			}
			if j+3 < len(tableNums[i]) && tableNums[i][j]*tableNums[i][j+1]*tableNums[i][j+2]*tableNums[i][j+3] > biggestHoriz {
				biggestHoriz = tableNums[i][j] * tableNums[i][j+1] * tableNums[i][j+2] * tableNums[i][j+3]
			}
			if i+3 < len(tableNums) && j+3 < len(tableNums[i]) && tableNums[i][j]*tableNums[i+1][j+1]*tableNums[i+2][j+2]*tableNums[i+3][j+3] > biggestDiag {
				biggestDiag = tableNums[i][j] * tableNums[i+1][j+1] * tableNums[i+2][j+2] * tableNums[i+3][j+3]
			}
			if i+3 < len(tableNums) && j+3 < len(tableNums[i]) && tableNums[i][j+3]*tableNums[i+1][j+2]*tableNums[i+2][j+1]*tableNums[i+3][j] > biggestRDiag {
				biggestRDiag = tableNums[i][j+3] * tableNums[i+1][j+2] * tableNums[i+2][j+1] * tableNums[i+3][j]
			}
		}
	}
	fmt.Println(biggestHoriz, biggestVert, biggestDiag, biggestRDiag)
}

func main() {
	problem11()
}
