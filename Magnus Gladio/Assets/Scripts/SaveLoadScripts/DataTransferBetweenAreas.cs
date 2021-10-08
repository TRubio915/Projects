using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public static class DataTransferBetweenAreas
{


    public static List<int> unlockedAbilities;
    public static List<int> equippedAbilitesMap;

    public static void UpdateUnlockedAbilities(List<int> abilities)
    {
        unlockedAbilities = abilities;
    }
    public static void UpdateEquippedAbilities(List<int> map)
    {
        equippedAbilitesMap = map;
    }



}
