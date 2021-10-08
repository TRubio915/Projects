using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public static class GameSettings
{
    public static float volumeLevel = 1;
    public static bool musicMuted = false;
    public static bool sfxMuted = false;

    public static bool isGamePaused = false;

    public static int currentScene = 1;
    public static int currentArea = 1;

    public static float[] playerPos;

    public static bool isNewGame = true;


    public static int saveSlot = 0;

    public static List<SaveData> saves = new List<SaveData>(); 


}
