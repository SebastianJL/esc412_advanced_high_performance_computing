#include <stdio.h>
#include <cstdint> // or (for C): #include <stdint.h>
#include <stdalign.h>


int main() {
    printf("sizeof(char)      = %lu\n", sizeof(char));
    printf("sizeof(short)     = %lu\n", sizeof(short));
    printf("sizeof(int)       = %lu\n", sizeof(int));
    printf("sizeof(long)      = %lu\n", sizeof(long));
    printf("sizeof(long long) = %lu\n", sizeof(long long));

    printf("sizeof(float)     = %lu\n", sizeof(float));
    printf("sizeof(double)    = %lu\n", sizeof(double));

    // This is C++ only
    printf("sizeof(bool)      = %lu\n", sizeof(bool));

    // From <cstdint> (C++) or <stdint.h> (C)
    printf("sizeof(int8_t)    = %lu\n", sizeof(std::int8_t));
    printf("sizeof(int16_t)   = %lu\n", sizeof(std::int16_t));
    printf("sizeof(int32_t)   = %lu\n", sizeof(std::int32_t));
    printf("sizeof(int64_t)   = %lu\n", sizeof(std::int64_t));

    struct s1 {
    	char a;
    	double b;
    	char c;
    	std::uint16_t d;
	} example;
    printf("sizeof(struct s1) = %lu\n", sizeof(struct s1));
    printf("sizeof(struct example.a) = %lu\n", sizeof(example.a));
    printf("sizeof(struct example.b) = %lu\n", sizeof(example.b));
    printf("sizeof(struct example.c) = %lu\n", sizeof(example.c));
    printf("sizeof(struct example.d) = %lu\n", sizeof(example.d));
    printf("alignof(struct s1) = %lu\n", alignof(s1));
    printf("alignof(struct example.a) = %lu\n", alignof(example.a));
    printf("alignof(struct example.b) = %lu\n", alignof(example.b));
    printf("alignof(struct example.c) = %lu\n", alignof(example.c));
    printf("alignof(struct example.d) = %lu\n", alignof(example.d));

    struct s2 {
    	char a;
    	char c;
    	std::uint16_t d;
    	double b;
	};
    printf("sizeof(struct s2) = %lu\n", sizeof(struct s2));
}
