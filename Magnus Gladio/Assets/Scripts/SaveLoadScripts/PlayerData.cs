using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class PlayerData : MonoBehaviour
{


    public Vector3 playerPosition;
    public int currentRegion;
    public int currentArea;
    public int playerCurr;
    public List<int> inventoryItems;
    public List<int> unlockedAbilitites; 
    public List<int> equippedAbilitiesMap;
    public int coinCount;

    public void Start()
    {
        inventoryItems = new List<int>();
    }
    public void UpdatePlayerPosition(Vector3 pos)
    {
        playerPosition = pos;
    }
    public void UpdateArea(int area)
    {
        currentArea = area;
    }
    public void updateRegion(int reg)
    {
        currentRegion = reg;
    }
    public void updatePlayerCurr(int curr)
    {
        playerCurr = curr;
    }

    public void updateInventory(List<int> items)
    {
        inventoryItems = items;
    }

    public void UpdateUnlockedAbilities(List<int> abilities)
    {
        unlockedAbilitites = abilities;
    }
    public void UpdateEquippedAbilities(List<int> map)
    {       
        equippedAbilitiesMap = map;
    }

    public void UpdateCoinCount(int count)
    {
        coinCount = count;
    }

}
