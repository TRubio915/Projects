using System.Collections.Generic;
using System.IO;
using System.Runtime.Serialization.Formatters.Binary;
using UnityEngine;
using UnityEngine.UIElements;

public static class SaveSystem
{
    public static void SaveGame(PlayerData player, int slot)
    {
        BinaryFormatter formatter = new BinaryFormatter();
        int currentSave = 0;     

        string path = Application.persistentDataPath + "/saveState" + slot + ".gp";
        FileStream fileStream = new FileStream(path, FileMode.Create);

        SaveData data = new SaveData(player);

        formatter.Serialize(fileStream, data);

        fileStream.Close();

    }

    public static SaveData LoadGame(int saveSlot)
    {
        string path = Application.persistentDataPath + "/saveState" + saveSlot +".gp";


        if(File.Exists(path))
        {
            BinaryFormatter formatter = new BinaryFormatter();

            FileStream fileStream = new FileStream(path, FileMode.Open);

            SaveData data = formatter.Deserialize(fileStream) as SaveData;
            fileStream.Close();

            return data;

        }
        else
        {
            //Debug.LogError("Save File Not found in " + path);
            return null;
        }   

    }

    public static List<string> getFileNames()
    {
        List<string> list = new List<string>();
       
        for(int i = 0; i < 5; i++)
        {
            string path = Application.persistentDataPath + "/saveState" + i + ".gp";

            if (File.Exists(path))
            {
                BinaryFormatter formatter = new BinaryFormatter();
                FileStream fileStream = new FileStream(path, FileMode.Open);
                SaveData data = formatter.Deserialize(fileStream) as SaveData;
                list.Add(data.saveTime.ToString());
                fileStream.Close();
            }
        }


        return list;
    }

}
