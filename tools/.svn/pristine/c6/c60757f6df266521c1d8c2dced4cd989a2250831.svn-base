import sys
import yaml

def get_value_from_yaml(era, key1, key2, key3):
    if era == '2018' : yaml_path =  "photon_prescales_2018.yaml"
    with open(yaml_path, 'r') as file:
        data = yaml.safe_load(file)
        try:
            return data[key1][key2][key3]
        except (KeyError, TypeError):
            return None

if __name__ == '__main__':
    if len(sys.argv) < 5:
        print("Usage: python script.py <yaml_path> <key1> <key2> <key3>")
        sys.exit(1)

    yaml_path = sys.argv[1]
    key1 = sys.argv[2]
    key2 = int(sys.argv[3])
    key3 = int(sys.argv[4])

    value = get_value_from_yaml("2018", yaml_path, key1, key2, key3)
    print(value)


