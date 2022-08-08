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

	case probNum == 12:
		//This is really slow
		var firstIndex int
		for i := 1; len(divisors(triangle(i))) <= 500; i++ {
			firstIndex = i
			fmt.Println(triangle(i), "has", len(divisors(triangle(i))), "divisors")
		}
		fmt.Println(triangle(firstIndex+1), "has", len(divisors(triangle(firstIndex+1))), "divisors")

	case probNum == 13:
		numStrings := []string{"37107287533902102798797998220837590246510135740250", "46376937677490009712648124896970078050417018260538", "74324986199524741059474233309513058123726617309629", "91942213363574161572522430563301811072406154908250", "23067588207539346171171980310421047513778063246676", "89261670696623633820136378418383684178734361726757", "28112879812849979408065481931592621691275889832738", "44274228917432520321923589422876796487670272189318", "47451445736001306439091167216856844588711603153276", "70386486105843025439939619828917593665686757934951", "62176457141856560629502157223196586755079324193331", "64906352462741904929101432445813822663347944758178", "92575867718337217661963751590579239728245598838407", "58203565325359399008402633568948830189458628227828", "80181199384826282014278194139940567587151170094390", "35398664372827112653829987240784473053190104293586", "86515506006295864861532075273371959191420517255829", "71693888707715466499115593487603532921714970056938", "54370070576826684624621495650076471787294438377604", "53282654108756828443191190634694037855217779295145", "36123272525000296071075082563815656710885258350721", "45876576172410976447339110607218265236877223636045", "17423706905851860660448207621209813287860733969412", "81142660418086830619328460811191061556940512689692", "51934325451728388641918047049293215058642563049483", "62467221648435076201727918039944693004732956340691", "15732444386908125794514089057706229429197107928209", "55037687525678773091862540744969844508330393682126", "18336384825330154686196124348767681297534375946515", "80386287592878490201521685554828717201219257766954", "78182833757993103614740356856449095527097864797581", "16726320100436897842553539920931837441497806860984", "48403098129077791799088218795327364475675590848030", "87086987551392711854517078544161852424320693150332", "59959406895756536782107074926966537676326235447210", "69793950679652694742597709739166693763042633987085", "41052684708299085211399427365734116182760315001271", "65378607361501080857009149939512557028198746004375", "35829035317434717326932123578154982629742552737307", "94953759765105305946966067683156574377167401875275", "88902802571733229619176668713819931811048770190271", "25267680276078003013678680992525463401061632866526", "36270218540497705585629946580636237993140746255962", "24074486908231174977792365466257246923322810917141", "91430288197103288597806669760892938638285025333403", "34413065578016127815921815005561868836468420090470", "23053081172816430487623791969842487255036638784583", "11487696932154902810424020138335124462181441773470", "63783299490636259666498587618221225225512486764533", "67720186971698544312419572409913959008952310058822", "95548255300263520781532296796249481641953868218774", "76085327132285723110424803456124867697064507995236", "37774242535411291684276865538926205024910326572967", "23701913275725675285653248258265463092207058596522", "29798860272258331913126375147341994889534765745501", "18495701454879288984856827726077713721403798879715", "38298203783031473527721580348144513491373226651381", "34829543829199918180278916522431027392251122869539", "40957953066405232632538044100059654939159879593635", "29746152185502371307642255121183693803580388584903", "41698116222072977186158236678424689157993532961922", "62467957194401269043877107275048102390895523597457", "23189706772547915061505504953922979530901129967519", "86188088225875314529584099251203829009407770775672", "11306739708304724483816533873502340845647058077308", "82959174767140363198008187129011875491310547126581", "97623331044818386269515456334926366572897563400500", "42846280183517070527831839425882145521227251250327", "55121603546981200581762165212827652751691296897789", "32238195734329339946437501907836945765883352399886", "75506164965184775180738168837861091527357929701337", "62177842752192623401942399639168044983993173312731", "32924185707147349566916674687634660915035914677504", "99518671430235219628894890102423325116913619626622", "73267460800591547471830798392868535206946944540724", "76841822524674417161514036427982273348055556214818", "97142617910342598647204516893989422179826088076852", "87783646182799346313767754307809363333018982642090", "10848802521674670883215120185883543223812876952786", "71329612474782464538636993009049310363619763878039", "62184073572399794223406235393808339651327408011116", "66627891981488087797941876876144230030984490851411", "60661826293682836764744779239180335110989069790714", "85786944089552990653640447425576083659976645795096", "66024396409905389607120198219976047599490197230297", "64913982680032973156037120041377903785566085089252", "16730939319872750275468906903707539413042652315011", "94809377245048795150954100921645863754710598436791", "78639167021187492431995700641917969777599028300699", "15368713711936614952811305876380278410754449733078", "40789923115535562561142322423255033685442488917353", "44889911501440648020369068063960672322193204149535", "41503128880339536053299340368006977710650566631954", "81234880673210146739058568557934581403627822703280", "82616570773948327592232845941706525094512325230608", "22918802058777319719839450180888072429661980811197", "77158542502016545090413245809786882778948721859617", "72107838435069186155435662884062257473692284509516", "20849603980134001723930671666823555245252804609722", "53503534226472524250874054075591789781264330331690"}
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

	case probNum == 31:
		britCoins := []int{1, 2, 5, 10, 20, 50, 100, 200}
		fmt.Println(len(coinCombos(britCoins, 200)))

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
		for i := 1; i < 10000000; i = i + 2 {
			bin, _ := strconv.Atoi(strconv.FormatInt(int64(i), 2))
			if isPalindrome(i) && isPalindrome(bin) {
				sum += i
				fmt.Println(i, bin)
			}
		}
		fmt.Println(sum)

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
			index := int(math.Pow(float64(10), float64(i)))
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

	default:
		fmt.Println("Haven't done this problem yet.")
	}
}
