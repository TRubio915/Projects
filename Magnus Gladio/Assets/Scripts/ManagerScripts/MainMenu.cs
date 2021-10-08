using System;
using System.Collections;
using System.Collections.Generic;
using TMPro;
using UnityEngine;
using UnityEngine.SceneManagement;
using UnityEngine.UI;

public class MainMenu : MonoBehaviour
{

    public Slider volumeSlider;
    public Toggle muteMusicToggle;
    public Toggle muteSFXToggle;
    public GameObject mainMenu;
    public Slider loadingBar;
    public TMP_Text progressText;

    [SerializeField]
    private TMP_Text[] saveslotsTexts;


    private void Start()
    {
        Cursor.visible = true;
        Time.timeScale = 1;
        GameSettings.isGamePaused = false;
        loadingBar.gameObject.SetActive(false);
        GameSettings.volumeLevel = volumeSlider.value;
    }

    public void LoadTutorial()
    {
        StartCoroutine(LoadAsynchronously(1));
        GameSettings.isNewGame = true;     
    }

    IEnumerator LoadAsynchronously(int sceneIndex)
    {
        AsyncOperation operation = SceneManager.LoadSceneAsync(sceneIndex);
        operation.allowSceneActivation = false;
        mainMenu.SetActive(false);
        loadingBar.gameObject.SetActive(true);

        while (loadingBar.value < 1)
        {
            float progress = Mathf.Clamp01(operation.progress / 0.9f);
            loadingBar.value = progress;
            progressText.text = (int) progress * 100f + "%";

            yield return null;
        }
        while (AudioManager.instance.soundsLoaded == false)
        {
            yield return null;
        }
        operation.allowSceneActivation = true;
    }

    public void DisplaySaves()
    {
        List<string> saveList = SaveSystem.getFileNames();

        if(saveList.Count > 0)
        {
            for(int i = 0; i < saveList.Count; i++)
            {
                saveslotsTexts[i].text = saveList[i];
            }           
        }
    }

    public void LoadGame(int slot)
    {
        SaveData saveData = SaveSystem.LoadGame(slot);
        if (saveData != null)
        {
            GameSettings.saveSlot = slot;
            GameSettings.isNewGame = false;
            GameSettings.currentArea = saveData.currentArea;
            GameSettings.currentScene = saveData.currentRegion;


            int areaToLoad = (saveData.currentRegion - 1) / 4 + 1;
            //GameSettings.playerPos = saveData.playerPosition;
            int sceneToLoad = (areaToLoad - 1) * 4 + 1;
            StartCoroutine(LoadAsynchronously(sceneToLoad));
        }
    }
    public void QuitGame()
    {
        Application.Quit();
    }
    public void VolumeChanged()
    {
        GameSettings.volumeLevel = volumeSlider.value;
        Sound[] sounds = FindObjectOfType<AudioManager>().sounds;

        foreach (Sound s in sounds)
        {
            s.source.volume = s.maxVolume * volumeSlider.value;
        }
    }

    public void MuteMusic()
    {
        GameSettings.musicMuted = muteMusicToggle.isOn;
        Sound[] sounds = FindObjectOfType<AudioManager>().sounds;
        Sound s = Array.Find(sounds, sound => sound.name == "Theme");
        if (s == null)
            return;
        s.source.mute = muteMusicToggle.isOn;
    }

    public void MuteSFX()
    {
        GameSettings.musicMuted = muteMusicToggle.isOn;
        Sound[] sounds = FindObjectOfType<AudioManager>().sounds;

        foreach(Sound s in sounds)
        {
            if (s == null)
                return;
            s.source.mute = muteSFXToggle.isOn;
        }

    }

}
