package main

import "fmt"

func main() {
	upperLimit := 1000
	var sum int
	for i := 0; i < upperLimit; i++ {
		if (i%3 == 0) || (i%5 == 0) {
			sum += i
		}
	}
	fmt.Println(sum)
}
