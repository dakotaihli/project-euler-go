package main

import (
	"fmt"
	"math"
	"math/big"
	"os"
	"sort"
	"strconv"
	"strings"
)

func problem(probNum int) {
	switch {
	case probNum == 1:
		upperLimit := 1000
		var sum int
		for i := 0; i < upperLimit; i++ {
			if (i%3 == 0) || (i%5 == 0) {
				sum += i
			}
		}
		fmt.Println(sum)

	case probNum == 2:
		upperLimit := 4000000
		var sum int
		for i := 0; fibonacci(i) < upperLimit; i++ {
			if i%3 == 1 {
				sum += fibonacci(i)
			}
		}
		fmt.Println(sum)

	case probNum == 3:
		N := 600851475143
		facs := primeFactorize(N)
		fmt.Println(facs[len(facs)-1])

	case probNum == 4:
		var found bool
		for i := 0; i < 1998 && !found; i++ {
			for j := 0; j <= i; j++ {
				if isPalindromeInt((999 - j) * (999 - (i - j))) {
					found = true
					fmt.Println(999-j, 999-(i-j), (999-j)*(999-(i-j)))
					break
				}
			}
			if found {
				break
			}
		}

	case probNum == 5:
		N := 20
		var firstN []int
		for i := 1; i <= N; i++ {
			firstN = append(firstN, i)
		}
		fmt.Println(lcm(firstN))

	case probNum == 6:
		N := 100
		fmt.Println(((N)*(N+1)/2)*((N)*(N+1)/2) - ((N)*(N+1)*(1+(2*N)))/6)

	case probNum == 7:
		N := 10001
		var primes []int
		for i := 2; len(primes) < N; i++ {
			if isPrime(i) {
				primes = append(primes, i)
			}
		}
		fmt.Println(primes[N-1])

	case probNum == 8:
		var numDigs []int
		dat, _ := os.ReadFile("p008_number.txt")
		bigNum := string(dat)
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

	case probNum == 9:
		var ag, bg, cg int
		for c := 334; c < 998; c++ {
			for b := 2; b+c < 999; b++ {
				if (b*b)+(1000-(b+c))*(1000-(b+c)) == c*c {
					ag, bg, cg = (1000 - (b + c)), b, c
				}
			}
		}
		fmt.Println(ag * bg * cg)

	case probNum == 10:
		n := 2000000
		primes := primesBelow(n)
		var sum int
		for _, p := range primes {
			sum += p
		}
		fmt.Println(sum)

	case probNum == 11:
		dat, _ := os.ReadFile("p011_table.txt")
		tableString := string(dat)
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

	case probNum == 12:
		//This is really slow
		var firstIndex int
		for i := 1; len(divisors(triangle(i))) <= 500; i++ {
			firstIndex = i
			fmt.Println(triangle(i), "has", len(divisors(triangle(i))), "divisors")
		}
		fmt.Println(triangle(firstIndex+1), "has", len(divisors(triangle(firstIndex+1))), "divisors")

	case probNum == 13:
		dat, _ := os.ReadFile("p013_numbers.txt")
		numStrings := strings.Split(string(dat), "\n")
		nums := make([]big.Int, len(numStrings))
		sum := new(big.Int)
		for i, x := range nums {
			x.SetString(numStrings[i], 10)
			sum.Add(sum, &x)
		}
		fmt.Println(sum.String())

	case probNum == 14:
		numberOfConcern := 1000000
		collatzLengths := map[int]int{
			1: 1,
		}
		for i := 2; i <= numberOfConcern; i++ {
			if _, ok := collatzLengths[i]; !ok {
				var collatzPath []int
				var defined bool
				for j := i; !defined; j = collatz(j) {
					collatzPath = append(collatzPath, j)
					_, defined = collatzLengths[j]
				}
				for j, N := 2, len(collatzPath); j <= N; j++ {
					collatzLengths[collatzPath[N-j]] = collatzLengths[collatzPath[N-1]] + j - 1
				}
			}
		}
		var indexLongest, lengthLongest int
		for n, L := range collatzLengths {
			if n < numberOfConcern && L > lengthLongest {
				indexLongest, lengthLongest = n, L
			}
		}
		fmt.Println(indexLongest, lengthLongest)

	case probNum == 15:
		choose := bigFact(20)
		choose.Mul(choose, choose)
		choose.Div(bigFact(40), choose)
		fmt.Println(choose.String())

	case probNum == 16:
		pow := new(big.Int)
		pow.Exp(big.NewInt(int64(2)), big.NewInt(int64(1000)), nil)
		fmt.Println(digitSum(pow.String()))

	case probNum == 17:
		var sum int
		for i := 1; i <= 1000; i++ {
			sum += len(strings.ReplaceAll(strings.ReplaceAll(numberName(i), " ", ""), "-", ""))
		}
		fmt.Println(sum)

	case probNum == 18 || probNum == 67:
		/* Consider a path down the triangle. If one wished to maximize the
		 * sum along the path, the choice at the last level is easy --- just
		 * the maximum of the left or right choices. Since at the penultimate
		 * level, one knows which option to take, we can replace that row
		 * by adding each cell with the maximum of its two children. The max
		 * sum of the resulting triangle is the same as the original.
		 * If n is the number of rows of the triangle, the operation of
		 * reducing in this way is O(n) (since n is also the number of items
		 * in the last level). Reducing down to a single row gives the maximum
		 * sum, and is O(n + (n-1) + ... + 2 + 1) = O(n^2).
		 */
		var filename string
		if probNum == 67 {
			filename = "p067_triangle.txt"
		} else {
			filename = "p018_triangle.txt"
		}
		dat, err := os.ReadFile(filename)
		if err != nil {
			panic(err)
		}
		var levels []string = strings.Split(string(dat), "\n")
		triStr := make([][]string, len(levels))
		tri := make([][]int, len(levels))
		for i := 0; i < len(levels); i++ {
			triStr[i] = strings.Split(levels[i], " ")
			tri[i] = make([]int, i+1)
			for j := 0; j <= i; j++ {
				tri[i][j], _ = strconv.Atoi(triStr[i][j])
			}
		}
		for len(tri) > 1 {
			tri, _ = trianglePathReduce(tri)
		}
		fmt.Println("The maximum path sum is", tri[0][0])

	case probNum == 19:
		var count int
		for y, m, d, dotW := 1900, 1, 1, 1; isDateBeforeOrSame(y, m, d, 2000, 12, 31); y, m, d = tomorrow(y, m, d) {
			if d == 1 && y >= 1901 && dotW == 0 {
				count++
			}
			dotW = (dotW + 1) % 7
		}
		fmt.Println(count)

	case probNum == 20:
		fmt.Println(digitSum(bigFact(100).String()))

	case probNum == 21:
		var sum int
		for a := 1; a < 10000; a++ {
			b := sumOfDivisors(a)
			if a != b && sumOfDivisors(b) == a {
				sum += a
			}
		}
		fmt.Println(sum)

	case probNum == 22:
		dat, err := os.ReadFile("p022_names.txt")
		if err != nil {
			panic(err)
		}
		var sum int
		names := strings.Split(strings.ReplaceAll(string(dat), "\"", ""), ",")
		sort.Strings(names)
		for i, s := range names {
			sum += (i + 1) * alphValue(s)
		}
		fmt.Println(sum)

	case probNum == 23:
		//This is a bit slow but not too bad
		var sum int
		for i := 2; i <= 28123; i++ {
			var canBe bool
			for j := 1; j < i; j++ {
				canBe = canBe || (isAbundant(j) && isAbundant(i-j))
			}
			if !canBe {
				sum += i
			}
		}
		fmt.Println(sum)

	case probNum == 24:
		//Count how many came before

	case probNum == 25:
		//Use the fact that F_n = floor(1/2 + (phi^n)/sqrt(5)), where phi is the golden ratio
		fmt.Println((999*math.Log(10) + math.Log(math.Sqrt(5))) / math.Log(math.Phi))

	case probNum == 26:
		/* Note that for any d \in \N, if 1/d = 0.(a1 a2 \ldots an), then we have 1/d =
		 * a1 a2 \ldots an / (10^n - 1). Since 2 and 5 are always relatively prime to
		 * 10^n - 1, multiplying each side by 2 or 5 does not cancel any factors on the
		 * right hand side, and hence the length of the cycle in the decimal of 2/d and
		 * 5/d stays the same. But since 2/d = 1/(d/2) and 5/d = 1/(d/5), we may divide
		 * out all factors of 2 and 5 from d to calculate the length of the cycles. This
		 * reduces the relevant calculations to those numbers relatively prime to 10.
		 *
		 * Then for all such numbers, from 1/d = a1 a2 \ldots an / (10^n - 1) it follows
		 * that n is simply the least such that d divides 10^n - 1.
		 */
		N := 1000
		var longestD, longestLength int
		length := make([]int, N)
		length[0], length[1] = 0, 0
		for d := 2; d < N; d++ {
			if d%2 == 0 {
				length[d] = length[d/2]
			} else if d%5 == 0 {
				length[d] = length[d/5]
			} else {
				length[d] = 1
				nines := big.NewInt(int64(9))
				bigD := big.NewInt(int64(d))
				mod := big.NewInt(int64(0))
				mod.Mod(nines, bigD)
				for int(mod.Int64()) != 0 {
					length[d]++
					nines.Add(nines, big.NewInt(int64(1)))
					nines.Mul(nines, big.NewInt(int64(10)))
					nines.Sub(nines, big.NewInt(int64(1)))
					mod.Mod(nines, bigD)
				}
			}
			if length[d] > longestLength {
				longestD, longestLength = d, length[d]
			}
		}
		fmt.Println("The reciprocal of", longestD, "has a decimal expansion with a cycle of length", longestLength)

	case probNum == 27:
		N := 1000
		maxConsecPrimes := 40
		maxPrimesCoeffA, maxPrimesCoeffB := 1, 41
		for b := -N; b <= N; b++ {
			if !isPrime(b) {
				// b is the value when n = 0, so for this to be prime, b must be prime
				continue
			}
			for a := 1 - N; a < N; a++ {
				consec := 0
				for n := 0; isPrime((n+a)*n + b); n++ {
					consec++
				}
				if consec > maxConsecPrimes {
					maxConsecPrimes = consec
					maxPrimesCoeffA, maxPrimesCoeffB = a, b
				}
			}
		}
		fmt.Println(maxPrimesCoeffA * maxPrimesCoeffB)

	case probNum == 28:
		// For k \geq 1, the (2k-1) by (2k-1) spiral can be thought of as layers radiating from 1. (The case k=1 is just the 1 by itself.)
		//For k \geq 2, the k-th layer is the numbers from (2k-3)^2 + 1 to (2k-1)^2, and the corners --- one of which is (2k-1)^2 --- are (2k-2) apart.
		//Working this out, the sum of the numbers in each layer is 16k^2 - 28k + 16. For a (2k-1) by (2k-1) spiral, the sum is the middle 1 plus the sum of each layer from 2 to k. Using sum formulas you can reduce this to a cubic polynomial in k.
		N := 1001
		k := (N + 1) / 2
		fmt.Println(-3 + 2*k - 6*k*k + (8*k*(1+2*k*k))/3)

	case probNum == 29:
		outs := make([]*big.Int, 0)
		for a := 2; a <= 100; a++ {
			for b := 2; b <= 100; b++ {
				bigA := big.NewInt(int64(a))
				bigB := big.NewInt(int64(b))
				c := big.NewInt(int64(0))
				c.Exp(bigA, bigB, nil)
				alreadyIn := false
				for _, x := range outs {
					//Somehow this bit makes it run really slow
					//TODO: Find more efficient way to tell if two bigInts are equal
					alreadyIn = alreadyIn || c.String() == x.String()
				}
				if !alreadyIn {
					outs = append(outs, c)
				}
			}
		}
		fmt.Println(len(outs))

	case probNum == 30:
		/* The largest that the sum of n-th powers can be is m*9^n, where
		 * m is the number of digits. Since 6*9^5 is only six digits, we
		 * know any of the desired numbers must have at most six digits.
		 * The problem also tells us to exclude the case 1 = 1^5, which
		 * would otherwise be the only one-digit solution. So we can bound
		 * our search by 10 <= N <= 6*9^5 = 354294
		 */
		uBound := 6 * intPow(9, 5)
		var sum int
		for N := 10; N <= uBound; N++ {
			digs := numToDigits(N)
			var digPowSum int
			for _, x := range digs {
				digPowSum += intPow(x, 5)
			}
			if digPowSum == N {
				sum += N
			}
		}
		fmt.Println("The sum of numbers that equal the fifth powers of their digits is", sum)

	case probNum == 31:
		britCoins := []int{1, 2, 5, 10, 20, 50, 100, 200}
		fmt.Println(len(coinCombos(britCoins, 200)))

	case probNum == 32:
		/* Without loss of generality, suppose the multiplier (the second
		 * number) is smaller than the multiplicand (the first number).
		 * The numbers of digits of the multiplier, multiplicand, and
		 * product must add to 9. This will provide the necessary bounds
		 * on our search
		 * If the multiplier (hence also the multiplicand) are three digits
		 * then the smallest that their product can be is 100*100 = 10000
		 * which makes a total of 3+3+5 = 11. Therefore, the multiplier
		 * must have at most two digits.
		 * If the multiplier has a single digit, and the multiplicand has
		 * three, then the largest the product can be is 9*999 = 8991 which
		 * makes a total of 3+1+4 = 8. Therefore, the multiplicand must
		 * have at least four digits whenever the multiplier has one.
		 * Indeed, it must have exactly four since the product must be
		 * larger, which means the only solutions must have 4+1+4 digits.
		 * If instead the multiplier has two digits, the multiplicand must
		 * have exactly three by similar reasoning.
		 * Furthermore, neither multiplier nor multiplcand may have a 1
		 * in the ones place, since then the ones digit of the other number
		 * will match the ones digit in the product, so the identity will
		 * not be pandigital.
		 * Similarly, no number may contain a 0.
		 */
		allDigits := []int{1, 2, 3, 4, 5, 6, 7, 8, 9}
		var solns [][]int
		for j := 2; j < 100; j++ {
			if j%10 == 0 || j%10 == 1 {
				continue
			}
			for i := 112; i*j < 10000; i++ {
				if i%10 == 0 || i%10 == 1 || (j < 10 && i < 1112) {
					continue
				}
				digits := append(append(numToDigits(i), numToDigits(j)...), numToDigits(i*j)...)
				sort.Ints(digits)
				if slicesEqual(digits, allDigits) {
					solns = append(solns, []int{i, j, i * j})
				}
			}
		}
		m := make(map[int]int)
		for _, s := range solns {
			m[s[2]]++
		}
		var sum int
		for x, _ := range m {
			sum += x
		}
		fmt.Println(sum)

	case probNum == 33:
		var sols [][]int
		for x := 1; x < 10; x++ {
			for y := x + 1; y < 10; y++ {
				for a := 1; a < 10; a++ {
					if (10*a+x)*y == (10*y+a)*x {
						sols = append(sols, []int{10*a + x, 10*y + a})
					}
					if (10*x+a)*y == (10*a+y)*x {
						sols = append(sols, []int{10*x + a, 10*a + y})
					}
				}
			}
		}
		fmt.Println(sols)
		num, denom := 1, 1
		for _, s := range sols {
			num *= s[0]
			denom *= s[1]
		}
		fmt.Println(denom / binGCD(num, denom))

	case probNum == 34:
		facs := []int{1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880}
		var sum int
		for n := 10; n <= 2540160; n++ {
			if (n < 100 && facs[n%10]+facs[(n/10)%10] == n) || (n >= 100 && n < 1000 && facs[n%10]+facs[(n/10)%10]+facs[(n/100)%10] == n) || (n >= 1000 && n < 10000 && facs[n%10]+facs[(n/10)%10]+facs[(n/100)%10]+facs[(n/1000)%10] == n) || (n >= 10000 && n < 100000 && facs[n%10]+facs[(n/10)%10]+facs[(n/100)%10]+facs[(n/1000)%10]+facs[(n/10000)%10] == n) || (n >= 100000 && n < 1000000 && facs[n%10]+facs[(n/10)%10]+facs[(n/100)%10]+facs[(n/1000)%10]+facs[(n/10000)%10]+facs[(n/100000)%10] == n) || (n >= 1000000 && facs[n%10]+facs[(n/10)%10]+facs[(n/100)%10]+facs[(n/1000)%10]+facs[(n/10000)%10]+facs[(n/100000)%10]+facs[(n/1000000)%10] == n) {
				sum += n
			}
		}
		fmt.Println(sum)

	case probNum == 35:
		circPrimes := []int{2, 3, 5, 7}
		for n := 11; n < 1000000; n++ {
			isCircPrime := isPrime(n)
			s := numToDigits(n)
			for i := 0; i < len(numToDigits(n)); i++ {
				s = append(s[1:], s[0])
				isCircPrime = isCircPrime && isPrime(digitsToNum(s))
			}
			if isCircPrime {
				circPrimes = append(circPrimes, n)
			}
		}
		fmt.Println("The", len(circPrimes), "circular primes are", circPrimes)

	case probNum == 36:
		var sum int
		N := 1000000
		for i := 1; i < N; i += 2 {
			//Only need to check odd i, since leading 0's don't count
			if !isPalindromeInt(i) {
				continue
			}
			binStr := strconv.FormatInt(int64(i), 2)
			if isPalindromeString(binStr) {
				sum += i
			}
		}
		fmt.Println(sum)

	case probNum == 37:
		var truncs []int
		var sum int
		//Use the given fact that there are 11 solutions to bound search
		for n := 11; len(truncs) < 11 && n < 1000000; n++ {
			if !isPrime(n) {
				continue
			}
			leftTrunc, rightTrunc := true, true
			for m := n / 10; m != 0; m /= 10 {
				rightTrunc = rightTrunc && isPrime(m)
			}
			for k, i := n%10, 1; intPow(10, i) <= n; k, i = n%intPow(10, i+1), i+1 {
				leftTrunc = leftTrunc && isPrime(k)
			}
			if leftTrunc && rightTrunc {
				truncs = append(truncs, n)
			}
		}
		fmt.Println(truncs)
		for _, n := range truncs {
			sum += n
		}
		fmt.Println(sum)

	case probNum == 38:
		/* Each value of n restricts the allowable sequence of digit-lengths,
		 * since the lengths of x, ..., nx have to add to 9. Also, note that
		 * the digits of x are the first digits of the concatenated product,
		 * which will help us restrict the search. Thus, we separate into
		 * cases for each value of n:
		 *
		 * If n >= 5, the only allowable digit lengths force x to be single-
		 * digit, but then the given example of 918273645 is the only
		 * solution with x = 9 (and hence the largest among n >=5 solutions).
		 *
		 * If n = 4, the only allowable digit lengths for x,2x,3x,4x must be
		 * 2,2,2,3 since 4x can have at most one more digit than x. But then
		 * if x is a 2-digit number such that 3x also has 2 digits, we have
		 * 10 <= x <= 33. In particular, the first digit of x (and hence the
		 * concatenated product also) must be 1, 2, or 3, making it smaller
		 * than the given solution of 918273645. Thus, we can ignore n = 4.
		 *
		 * If n = 3, by similar reasoning we must have x, 2x, and 3x are all
		 * 3-digit numbers. Then we must have 100 <= x <= 333, so any
		 * solution would also be less than the one given in the problem.
		 *
		 * Finally, if n = 2, then x has 4 digits and 2x has 5. Thus, 5000
		 * <= x <= 9999. But since we are looking for a solution larger
		 * than 918273645, and the result must be pandigital, we can further
		 * restrict the search to 9182 <= x <= 9876. Moreover, the pandigital
		 * product must be of the form 100002 * x.
		 */
		var sols []int
		for x := 9182; x <= 9876; x++ {
			if areNumsPerms(123456789, 100002*x) {
				sols = append(sols, 100002*x)
			}
		}
		fmt.Println(sols)

	case probNum == 39:
		trips := make([]int, 1001)
		for m := 2; m <= 22; m++ {
			for n := 1; n < m; n++ {
				if (m+n)%2 == 1 && binGCD(m, n) == 1 {
					for k := 1; k*m*(m+n) <= 500; k++ {
						trips[2*k*m*(m+n)]++
					}
				}
			}
		}
		var maxTriples, maxIndex int
		for i, s := range trips {
			if s > maxTriples {
				maxIndex, maxTriples = i, s
			}
		}
		fmt.Println("There are", maxTriples, "Pythagorean triples that sum to", maxIndex)

	case probNum == 40:
		champ := "."
		for i := 1; len(champ) <= 1000001; i++ {
			champ = champ + strconv.Itoa(i)
		}
		prod := 1
		for i := 0; i <= 6; i++ {
			index := intPow(10, i)
			dig, _ := strconv.Atoi(champ[index : index+1])
			prod *= dig
		}
		fmt.Println(prod)

	case probNum == 41:
		var pandigitalPrimes []int
		firstDigits := []int{1, 2, 3, 4, 5, 6, 7}
		for n := 1; n < 10000000; n++ {
			nDigits := numToDigits(n)
			sort.Ints(nDigits)
			if slicesEqual(nDigits, firstDigits[:len(nDigits)]) && isPrime(n) {
				pandigitalPrimes = append(pandigitalPrimes, n)
			}
		}
		fmt.Println(pandigitalPrimes)

	case probNum == 42:
		dat, err := os.ReadFile("p042_words.txt")
		if err != nil {
			panic(err)
		}
		words := strings.Split(strings.ReplaceAll(string(dat), "\"", ""), ",")
		triWords := 0
		for _, s := range words {
			for i := 1; triangle(i) <= alphValue(s); i++ {
				if triangle(i) == alphValue(s) {
					triWords++
				}
			}
		}
		fmt.Println(triWords)

	case probNum == 43:
		// The indents make this code miserable to read. Is there a better way?
		digs := []int{0, 1, 2, 3, 4, 5, 6, 7, 8, 9}
		var sols []int
		for i1, d1 := range digs {
			//if d1 == 0 {
			if false {
				continue
			} else {
				digsLess1 := removeSliceElt(digs, i1)
				for i2, d2 := range digsLess1 {
					digsLess2 := removeSliceElt(digsLess1, i2)
					for i3, d3 := range digsLess2 {
						digsLess3 := removeSliceElt(digsLess2, i3)
						for i4, d4 := range digsLess3 {
							if d4%2 != 0 {
								continue
							} else {
								digsLess4 := removeSliceElt(digsLess3, i4)
								for i5, d5 := range digsLess4 {
									if (d3+d4+d5)%3 != 0 {
										continue
									} else {
										digsLess5 := removeSliceElt(digsLess4, i5)
										for i6, d6 := range digsLess5 {
											if d6%5 != 0 {
												continue
											} else {
												digsLess6 := removeSliceElt(digsLess5, i6)
												for i7, d7 := range digsLess6 {
													if ((100*d5)+(10*d6)+d7)%7 != 0 {
														continue
													} else {
														digsLess7 := removeSliceElt(digsLess6, i7)
														for i8, d8 := range digsLess7 {
															if d6+d8-d7%11 != 0 {
																continue
															} else {
																digsLess8 := removeSliceElt(digsLess7, i8)
																for i9, d9 := range digsLess8 {
																	if ((100*d7)+(10*d8)+d9)%13 != 0 {
																		continue
																	} else {
																		d10 := digsLess8[1-i9] // There should only be one more now
																		if ((100*d8)+(10*d9)+d10)%17 == 0 {
																			sols = append(sols, digitsToNum([]int{d1, d2, d3, d4, d5, d6, d7, d8, d9, d10}))
																		}
																	}
																}
															}
														}
													}
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		fmt.Println(sols)
		var sum int
		for _, x := range sols {
			sum += x
		}
		fmt.Println("Sum is", sum)

	case probNum == 44:
		/* Suppose we have found k,j such that P_k - P_j and P_k + P_j are both
		 * pentagonal. Then the desired minimum of all such |P_k' - P_j'| is
		 * bounded above by |P_k - P_j|. This lets us put an upper bound on the
		 * values of k' and j' that minimize the desired quantity: the smallest
		 * distance from P_m to any other pentagonal distance is P_(m-1), and
		 * we can simply calculate that |P_m - P_(m-1)| = 3m-2. Thus, if we
		 * choose m such that |P_k - P_j| < 3m-2, we can then simply brute force
		 * all k',j' <= m to find the minimizing pair.
		 */
		pents := []int{1}
		var maxIndex int
		var minDist int
		for k, j := 2, 1; maxIndex == 0; k, j = pairIterate(k, j) {
			if k == 0 || j == 0 || k <= j {
				continue
			}
			for pents[len(pents)-1] < pentagonal(k)+pentagonal(j) {
				pents = append(pents, pentagonal(len(pents)+1))
			}
			if sortedIntSliceContains(pents, pents[k-1]+pents[j-1]) && sortedIntSliceContains(pents, pents[k-1]-pents[j-1]) {
				maxIndex = ((k - j) * (k + j)) / 2
				minDist = pents[k-1] - pents[j-1]
			}
		}
		for pents[len(pents)-1] < 2*pentagonal(maxIndex) {
			pents = append(pents, pentagonal(len(pents)+1))
		}
		fmt.Println(maxIndex, len(pents), minDist)
		// The following code finds the minimizing pair by brute force.
		// It's not good since the max index is about 1.8 million.
		// As it turns out, the pair found by the above code is already the minimizer.
		// I have no proof of this other than PE accepts the answer.
		// TODO: find a fast program that verifies this solution is correct.
		/*		for k := 2; k <= maxIndex; k++ {
					for j := 1; j < k; j++ {
						if sortedIntSliceContains(pents, pents[k-1]+pents[j-1]) {
							if sortedIntSliceContains(pents, pents[k-1]-pents[j-1]) {
								minDist = min(pents[k-1]-pents[j-1], minDist)
							}
						}
					}
				}
				fmt.Println("Among pairs of pentagonal numbers whose sum and difference are also pentagonal,", minDist, "is the smallest distance between them")
		*/

	case probNum == 45:
		/* Observe that H_n = T_(2n-1) for all n. Thus, it suffices to find the
		 * next hexagonal number which is also pentagonal.
		 * Calculation shows that any solution to T_n = P_m is given by m =
		 * (1+sqrt(48n^2 - 24n + 1))/6, which is bounded below by the easier
		 * sqrt((4/3)n^2 - (2/3)n). Thus we only need to check a couple m's
		 * per n.
		 */
		for n := 144; true; n++ {
			var isPentag bool
			for m := int(math.Floor((math.Sqrt((4.0/3.0)*math.Pow(float64(n), 2) - (2.0/3.0)*float64(n))))); pentagonal(m) <= hexagonal(n); m++ {
				isPentag = pentagonal(m) == hexagonal(n)
			}
			if isPentag {
				fmt.Println("Next hexagonal which is pentagonal is H_", n, "=", hexagonal(n))
				break
			}
		}

	case probNum == 46:
		for n := 9; true; n += 2 {
			isGold := false
			if !isPrime(n) {
				for i := 1; n > 2*i*i; i++ {
					if isPrime(n - (2 * i * i)) {
						isGold = true
						break
					}
				}
			}
			if !(isGold || isPrime(n)) {
				fmt.Println("Smallest counterexample is", n)
				break
			}
		}

	case probNum == 48:
		sum := new(big.Int)
		pow := new(big.Int)
		tens := big.NewInt(int64(10000000000))
		for i := 1; i <= 1000; i++ {
			pow.Exp(big.NewInt(int64(i)), big.NewInt(int64(i)), tens)
			sum.Add(sum, pow)
		}
		sum.Mod(sum, tens)
		fmt.Println(sum.String())

	case probNum == 49:
		var sols [][]int
		for N := 1000; N < 10000 && len(sols) < 2; N++ {
			if isPrime(N) {
				for d := 2; N+(2*d) < 10000; d += 2 {
					if isPrime(N+d) && isPrime(N+(2*d)) && areNumsPerms(N, N+d) && areNumsPerms(N, N+(2*d)) {
						sols = append(sols, []int{N, d})
					}
				}
			}
		}
		fmt.Println(sols)

	case probNum == 50:
		N := 1000000
		primesBelow1Mil := primesBelow(N)
		/* Let M be the largest number such that the first M-many primes (starting at 2)
		 * add up to at most N. This is a decent upper bound for the desired answer to the
		 * problem. Indeed, any longer sequence must add up to at least the sum of the
		 * first M+1 primes, which is strictly above N.
		 *
		 * longestPossible is this M.
		 */
		longestPossible := 0
		leastSum := 0
		for i := 0; leastSum <= N; i++ {
			leastSum += primesBelow1Mil[i]
			longestPossible++
		}
		var longestSeq []int
		var biggestPrime int
		var found bool
		for l := longestPossible; l >= 21; l-- {
			for i := 0; i+l <= len(primesBelow1Mil); i++ {
				if l*primesBelow1Mil[i] > N {
					continue
				}
				sum := 0
				for _, p := range primesBelow1Mil[i : i+l] {
					sum += p
				}
				if sum < N && isPrime(sum) {
					found = true
					longestSeq = primesBelow1Mil[i : i+l]
					biggestPrime = sum
					break
				} else if !isPrime(sum) {
					fmt.Println("The sum of", l, "consecutive primes starting at", primesBelow1Mil[i], "is", sum, " which is not prime")
				}
			}
			if found {
				break
			}
		}
		fmt.Println(len(longestSeq), "consecutive primes add to", biggestPrime)

	case probNum == 52:
		for i := 1; true; i++ {
			isGood := true
			numDigs := numToDigits(i)
			sort.Ints(numDigs)
			for n := 2; n <= 6 && isGood; n++ {
				mulDigs := numToDigits(n * i)
				sort.Ints(mulDigs)
				isGood = isGood && slicesEqual(numDigs, mulDigs)
			}
			if isGood {
				fmt.Println(i, "is good")
				break
			}
		}
		//TODO: The answer to 52 is 142857, which is 999999/7. Investigate?

	case probNum == 53:
		N, M := 100, 1000000
		//There are triangle(n+1)-1 possible values of nCr for n \geq 1
		//So we count the ones below 1 million and subtract
		var count int
		for n := 1; n <= N; n++ {
			for r := 0; 2*r < n && nCr(n, r) < M; r++ {
				count += 2 //Counting both nCr and nC(n-r)
			}
			if n < 23 && n%2 == 0 {
				count++ //Counting nC(n/2) for even n
			}
		}
		fmt.Println("There are", triangle(N+1)-1-count, "values above", M)

	case probNum == 63:
		var nums [][]int
		for n := 1; n <= 22; n++ {
			for i := 1; int(float64(n)*math.Log10(float64(i))) <= n-1; i++ {
				if int(float64(n)*math.Log10(float64(i))) == n-1 {
					nums = append(nums, []int{i, n})
				}
			}
		}
		fmt.Println(len(nums))

	case probNum == 293:
		N := 1000000000
		var admissibles []int
		for n := 2; n < N; n += 2 {
			if isAdmissible(n) {
				admissibles = append(admissibles, n)
			}
		}
		pseuForts := make(map[int][]int)
		for _, n := range admissibles {
			pF := primeAfter(n+1) - n
			pseuForts[pF] = append(pseuForts[pF], n)
		}
		var sum int
		for pF, _ := range pseuForts {
			sum += pF
		}
		fmt.Println("The sum of distinct pseudo-Fortunate numbers is", sum)

	default:
		fmt.Println("Haven't done this problem yet.")
	}
}
