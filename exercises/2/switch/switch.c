#include <stdio.h>

void name_digit(int i){
	switch(i){
		case 0:
			printf("Zero\n");
			break;

		case 1:
			printf("One\n");
			break;

		case 2:
			printf("Two\n");
			break;

		case 3:
			printf("Three\n");
			break;

		case 4:
			printf("Four\n");
			break;

		case 5:
			printf("Five\n");
			break;

		case 6:
			printf("Six\n");
			break;

		case 7:
			printf("Seven\n");
			break;

		case 8:
			printf("Eight\n");
			break;

		case 9:
			printf("Nine\n");
			break;

		default:
		printf("Not a digit\n");

	}
}

int main()
{
	name_digit(0);
	name_digit(1);
	name_digit(2);
	name_digit(3);
	name_digit(4);
	name_digit(5);
	name_digit(6);
	name_digit(7);
	name_digit(8);
	name_digit(9);
	return 0;
}