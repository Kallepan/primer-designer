// Group an Map with an array of objects by a property
export const groupBy = <T,>(array: T[], property: string) => {
    const map = new Map<string, T[]>();
    return array.reduce((_: any, currentValue: any) => {
      // Get the value of the property we want to group by
      const key = currentValue[property];
      
      // Get the current list of values for this key
      const list = map.get(key);
      if (list === undefined) {
        // If there is no list, create one
        map.set(key, [currentValue]);
      } else {
        // If there is a list, add the current value to it
        list.push(currentValue);
      }
      return map;
    }, new Map<string, T[]>());
  };
  