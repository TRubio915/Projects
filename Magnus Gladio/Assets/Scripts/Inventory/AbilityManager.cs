using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class AbilityManager : MonoBehaviour
{


    public GameObject[] UnlockedAbilities;
    public GameObject[] UnlockScreens;    
    public GameObject[] EquippedAbilities;
    public Dictionary<GameObject, int> EquippedAbilitiesDictionary;

    public int MaxEquippedAbilities;

    private GameManager manager; 
    private AbilityManager instance;

    private void Awake()
    {

        if (instance == null)
            instance = this;
        else
        {
            Destroy(gameObject);
            return;
        }
      //  DontDestroyOnLoad(this.gameObject);
        Setup();
    }

    public void Setup()
    {
        EquippedAbilitiesDictionary = new Dictionary<GameObject, int>();

        for (int i = 0; i < MaxEquippedAbilities; i++)
        {
            EquippedAbilities[i].GetComponent<Image>().enabled = false;
            EquippedAbilities[i].GetComponent<Button>().interactable = false;
            EquippedAbilitiesDictionary.Add(EquippedAbilities[i], -1);
        }

        manager = GetComponent<GameManager>();
    }

    public void equipAbility(int unlockedAbilityIndex)
    {

        if(unlockedAbilityIndex < UnlockedAbilities.Length)
        {
            for(int i = 0; i < MaxEquippedAbilities; i++)
            {
                if(EquippedAbilities[i].GetComponent<Image>().enabled == false)
                {
                    //EquippedAbilities[i] = UnlockedAbilities[unlockedAbilityIndex];
                    EquippedAbilities[i].GetComponent<Image>().enabled = true;
                    EquippedAbilities[i].GetComponent<Image>().sprite = UnlockedAbilities[unlockedAbilityIndex].GetComponent<Image>().sprite;
                    EquippedAbilities[i].GetComponent<Button>().interactable = true;
                    EquippedAbilitiesDictionary[EquippedAbilities[i]] = unlockedAbilityIndex;
                    UnlockedAbilities[unlockedAbilityIndex].GetComponent<Button>().interactable = false;

                    manager.EnableAbilityInGame(UnlockedAbilities[unlockedAbilityIndex].name);


                    break;
                }                
            }
        }

    }

    public void unEquipAbility(int equippedAbilityIndex)
    {
        GameObject equippedAbility = EquippedAbilities[equippedAbilityIndex];
        EquippedAbilities[equippedAbilityIndex].GetComponent<Image>().enabled = false;
        EquippedAbilities[equippedAbilityIndex].GetComponent<Button>().interactable = false;
        
        int unlockedAbilityIndex = EquippedAbilitiesDictionary[equippedAbility];
        UnlockedAbilities[unlockedAbilityIndex].GetComponent<Button>().interactable = true;

        manager.DisableAbilityInGame(UnlockedAbilities[unlockedAbilityIndex].name);
    }
  

    public void UnlockAbility(string abilityName)
    {
        for(int i = 0; i < UnlockedAbilities.Length; i++)
        {
            if(UnlockedAbilities[i].name == abilityName)
            {
                UnlockedAbilities[i].GetComponent<Image>().enabled = true;
                UnlockedAbilities[i].GetComponent<Button>().interactable = true;
                StartCoroutine(ShowUnlockScreen(i));
                break;
            }
        }
    }

    private IEnumerator ShowUnlockScreen(int index)
    {
        if (GameSettings.isNewGame)
        {
            yield return new WaitForSeconds(0.1f);
            Cursor.visible = true;
            GameObject.FindGameObjectWithTag("Player").GetComponent<Animator>().updateMode = AnimatorUpdateMode.Normal;

            Time.timeScale = 0;
            manager.isGamePaused = true;
            GameSettings.isGamePaused = true;
            UnlockScreens[index].SetActive(true);
        }
    }

    public void closeUnlockScreen()
    {

        Time.timeScale = 1;
        GameObject.FindGameObjectWithTag("Player").GetComponent<Animator>().updateMode = AnimatorUpdateMode.UnscaledTime;
        Cursor.visible = false;
        manager.isGamePaused = false;
        GameSettings.isGamePaused = false;
    }

    public void ReUnlockAbility(string abilityName)
    {
        for (int i = 0; i < UnlockedAbilities.Length; i++)
        {
            if (UnlockedAbilities[i].name == abilityName)
            {
                UnlockedAbilities[i].GetComponent<Image>().enabled = true;
                UnlockedAbilities[i].GetComponent<Button>().interactable = true;
                //StartCoroutine(ShowUnlockScreen(i));
                break;
            }
        }
    }

}
