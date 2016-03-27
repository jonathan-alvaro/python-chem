def loadMassList(filepath="massList.txt"):
    """
    Returns a dictionary that maps each element to its molar mass

    Data is loaded from filepath
    """
    input_file = open(filepath)
    
    massListDict = {}
    
    for line in input_file:

        line = line.split()
        if line[0] == 'Element' or line[0][0] == '-':
            continue

        elem_name = line[0]
        molar_mass = float(line[1])

        massListDict[elem_name] = molar_mass

    return massListDict

def polyIon(formula, amount=1):
    """
    Maps each element in a polyIon to its quantity in said ion. Then, multiply
    by amount (consider polyIon part of a molecule in which there exists AMOUNT
    molecules of polyIon

    Returns a dictionary
    """
    element_dict = readFormula(formula)

    for key in element_dict:
        element_dict[key] *= amount
    
    return element_dict

def readFormula(formula):
    """
    Reads in a molecular formula and returns a dictionary mapping each element
    in the formula to its quantity in every molecule of the formula given

    Assumes formula is a string and formula must be written correctly
    (i.e. element symbols are written with uppercase and a trailing lowercase,
    if any, etc.)
    """

    assert len(formula) > 0, "Formula is an empty string"

    current_element = formula[0]
    current_amount = 1
    formula_dict = {}

    index = 1

    #If first letter is an opening parentheses, then go back to index 0 to make
    #sure it's processed correctly
    if current_element == '(':
        index = 0
    
    while index < len(formula):
        
        letter = formula[index]
        
        if letter == '(':
            print "In polyion, index: " + str(index)
            poly_ion = {}
            ion_formula = ''

            #Increment index so it evaluates the next letter after '(', not the
            #parentheses itself
            index += 1

            while formula[index] != ')':

                #Asserts that end of formula has not been reached
                assert index < len(formula), "Missing closing parentheses!"

                ion_formula += formula[index]

                index += 1

            #After exiting, check if next character in formula is a number, if
            #yes, process it for polyIon functioin, if not leave it alone
            if formula[index + 1].isdigit():
                poly_ion = polyIon(ion_formula, int(formula[index+1]))
                index += 1
            
            #Then, after each element in polyIon has been accounted for, add it
            #to the whole formula's count
            for key in poly_ion:
                if key in formula_dict:
                    formula_dict[key] += poly_ion[key]

                else:
                    formula_dict[key] = poly_ion[key]
    
        #If letter is lowercase simply append to current element name
        elif letter.islower():

            current_element += letter

        #If letter is another uppercase, then that means there is only 1 of
        #current_element
        elif letter.isupper():
            
            #If there was no previous occurence of this element, add it to the
            #dictionary
            if current_element not in formula_dict:
                formula_dict[current_element] = current_amount

            #Simply increase the quantity if element had been accounted for
            else:
                formula_dict[current_element] += current_amount

            current_element = letter

        #If it is a number, then change the current_amount and add it to the
        #dictionary. Then, reset the variables
        elif letter.isdigit():

            current_amount = int(letter)

            #If the digit in the formula is longer than 1 digit and it is not
            #the last character in the formula
            while index < (len(formula) - 1) and formula[index+1].isdigit():
                current_amount *= 10
                current_amount += int(formula[index+1])

                index += 1

            if current_element not in formula_dict:
                formula_dict[current_element] = current_amount

            else:
                formula_dict[current_element] += current_amount
                
            current_element = ''
            current_amount = 1

        else:
            raise ValueError, "formula contains invalid character"
        
        index += 1

    #Makes sure final element in formula is added to dictionary
    if (current_element not in formula_dict) and current_element.isalpha():
        formula_dict[current_element] = current_amount

    elif current_element.isalpha():
        formula_dict[current_element] += current_amount

    #Remove all invalid keys in the formula, i.e. parentheses and empty strings
    for key in formula_dict.keys():
        if not key.isalpha():
            formula_dict.pop(key)
        
    return formula_dict

def calculateMass(formula):
    """
    Calculates the molar mass of a chemical compound using its formula

    Assumes formula is a molecular formula written following the standard IUPAC
    rules

    Returns a float
    """
    element_count = readFormula(formula)

    massList = loadMassList()

    molar_mass = 0
    
    for element in element_count:
        molar_mass += massList[element] * element_count[element]

    return round(molar_mass, 2)

def gcd(num1, num2):
    """
    Returns the greatest common divisor of two numbers

    Assumes num1 and num2 are numeric datas

    Uses Binary GCD Algorithm, for more info:
    https://en.wikipedia.org/wiki/Binary_GCD_algorithm

    Returns a float of a whole number
    """

    if num1 == 0 or num2 == 0:
        return float(max(num1, num2))

    elif num1 % 2 == 0:

        if num2 % 2 == 0:
            return 2 * gcd(num1/2, num2/2)

        else:
            return gcd(num1/2, num2)

    elif num2 % 2 == 0:

        if num1 % 2 == 0:
            return 2 * gcd(num1/2, num2/2)

        else:
            return gcd(num1, num2/2)

    else:
        bigger = max(num1, num2)
        smaller = min(num1, num2)

        return gcd((bigger-smaller)/2, smaller)
