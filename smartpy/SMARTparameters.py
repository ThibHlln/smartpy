from csv import DictReader


def get_parameters_from_file(file_location):
    my_dict_par = dict()
    try:
        with open(file_location, 'rb') as my_file:
            my_reader = DictReader(my_file)
            for row in my_reader:
                my_dict_par[row['PAR_NAME']] = float(row['PAR_VALUE'])
    except KeyError:
        raise Exception("There is 'PAR_NAME' or 'PAR_VALUE' column in {}.".format(file_location))
    except ValueError:
        raise Exception("There is at least one incorrect parameter value in {}.".format(file_location))
    except IOError:
        raise Exception("There is no parameters file at {}.".format(file_location))

    return my_dict_par
