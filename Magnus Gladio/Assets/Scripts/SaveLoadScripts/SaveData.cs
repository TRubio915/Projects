using System.Collections;
using System.Collections.Generic;
using UnityEditor;
using UnityEngine;
using System;

[System.Serializable]
public class SaveData
{

    public float playerHealth;
    public float playerStamina;
    public float[] playerPosition;
    public int currentRegion;
    public int currentArea;
    public List<int> inventoryItems;
    public List<int> unlockedAbilities;
    public List<int> equippedAbilititesMap;
    public int coinCount;


    public DateTime saveTime;
    public float gameVolumeLevel;


    public SaveData(PlayerData player)
    {
        saveTime = System.DateTime.Now;
        playerHealth = PlayerStats.PlayerHealth;
        playerStamina = PlayerStats.PlayerStamina;

        playerPosition = new float[3];

        playerPosition[0] = player.playerPosition.x;
        playerPosition[1] = player.playerPosition.y;
        playerPosition[2] = player.playerPosition.z;
        currentRegion = player.currentRegion;
        currentArea = player.currentArea;

        /*
        inventoryItems = new List<int>();
        inventoryItems = player.inventoryItems;
        */

        unlockedAbilities = new List<int>();
        unlockedAbilities = player.unlockedAbilitites;

        equippedAbilititesMap = new List<int>();
        equippedAbilititesMap = player.equippedAbilitiesMap;


        coinCount = player.coinCount;

        gameVolumeLevel = GameSettings.volumeLevel;
       
    }

}
