import random
import string

def get_random_string(length):
    # Random string with the combination of lower and upper case
    letters = string.ascii_letters
    result_str = ''.join(random.choice(letters) for i in range(length))
    f = open("random_strings.txt", "a")
    f.write(result_str)
    f.write("\n");
    f.close()
    #print("Random string is: ", result_str)

for i in range (1, 100000):
    get_random_string(1000)