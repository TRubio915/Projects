
using JetBrains.Annotations;
using System.Collections;
using System.Collections.Generic;
using TMPro;
using UnityEngine;
using UnityEngine.SceneManagement;
using UnityEngine.UI;

public class GameManager : MonoBehaviour
{

    public GameObject ingameMenu;
    public GameObject saveMenu;
    public GameObject Inventory;
    public GameObject camObject;
    public Upgrademanager upgradeManager;
    public Image throwableIndicator;
    public Slider healthBarSlider;
    public Slider staminaBarSlider, staminaBarWhiteSlider;
    public HealthBar healthBarScript;
    public Slider volumeSlider;

    private PlayerController playerController;
    private RangedAttackWithTeleportation r_Controller;


    public int coinsCollected = 0;
    public int maxCoinsAllowed;
    public TMP_Text coinCount;


    [SerializeField]
    private Animator sceneTransition;

    public float 
        staminaCostToTeleport,
        staminaCostToRun,
        staminaRegenSpeed;

    private int currentSceneNumber = 1, currentArea = 1, nextSceneNumber, previousSceneNumber, regionPart = 0;

    public bool isGameStart = false;
    private bool isInventoryOpen;

    public bool isGamePaused = false;

    
    [SerializeField]
    private GameObject healthPotion;

    [SerializeField]
    private List<TMP_Text> saveslotsTexts;

    public Image healthCircleBar;


    [Space(10)]
    [Header("Tutorial")]

    public bool isFirstGame = false;
    public bool tutorialBegin = true;
    public bool isFirstTimeInArea4 = true;

    public GameObject TutorialUIElements;

    public GameObject KunaiCoolDownTutorial;

    private GameManager instance;

    [SerializeField]
    private GameObject mode2;
    public bool isKunaiModeUnlocked = false;
    private bool isKunaiModeTutFinished = false;
    [SerializeField]
    private GameObject kunaiModeTutObj;

    private void Awake()
    {

        if (instance == null)
            instance = this;
        else
        {
            Destroy(gameObject);
            return;
        }
        //DontDestroyOnLoad(this.gameObject);

    }

    // Start is called before the first frame update
    void Start()
    {
        Time.timeScale = 1;
        currentSceneNumber = SceneManager.GetActiveScene().buildIndex;
        currentArea = (currentSceneNumber - 1) / 4 + 1;
        StartCoroutine(InitiateLevel());
    }

    IEnumerator InitiateLevel()
    {

        if (!GameSettings.isNewGame && (GameSettings.currentScene - 1 ) % 4 > 0 )
        {
            AsyncOperation operation = SceneManager.LoadSceneAsync(GameSettings.currentScene, LoadSceneMode.Additive);
            while (!operation.isDone)
            {
                yield return null;
            }
        }

        GameObject sceneTransitionSprite = GameObject.Find("SceneTransitionSprite");

        sceneTransitionSprite.GetComponent<SpriteRenderer>().enabled = true;
        sceneTransition = sceneTransitionSprite.GetComponent<Animator>();


        playerController = GameObject.FindGameObjectWithTag("Player").GetComponent<PlayerController>();
        r_Controller = playerController.GetComponent<RangedAttackWithTeleportation>();
        PlayerStats.UpdateMaxPlayerHealth(100);
        PlayerStats.UpdateMaxPlayerStamina(100);
        PlayerStats.UpdatePlayerHealth(100);
        PlayerStats.UpdatePlayerStamina(100);
        healthBarScript = GameObject.FindObjectOfType<HealthBar>();
        isInventoryOpen = false;
        isGameStart = true;


        if (SaveSystem.getFileNames().Count <= 0)
        {
            isFirstGame = true;
            StartCoroutine(BeginUITutorial());
        }
        else
        {            
            if (!GameSettings.isNewGame)
            {
                LoadGame(GameSettings.saveSlot);                    
            }
        }

        AbilityManager manager = GetComponent<AbilityManager>();

        manager.Setup();

        if (DataTransferBetweenAreas.unlockedAbilities != null)
        {
            for (int i = 0; i < DataTransferBetweenAreas.unlockedAbilities.Count; i++)
            {
                if (DataTransferBetweenAreas.unlockedAbilities[i] != -1)
                {
                    int index = DataTransferBetweenAreas.unlockedAbilities[i];
                    GameObject obj = manager.UnlockedAbilities[index];
                    manager.ReUnlockAbility(obj.name);
                }
            }

            for (int i = 0; i < DataTransferBetweenAreas.equippedAbilitesMap.Count; i++)
            {

                int mapIndex = DataTransferBetweenAreas.equippedAbilitesMap[i];
                if (mapIndex != -1)
                {
                    manager.equipAbility(mapIndex);
                }
            }

        }

        Sound[] sounds = FindObjectOfType<AudioManager>().sounds;
        foreach (Sound s in sounds)
        {
            s.source = gameObject.AddComponent<AudioSource>();
            s.source.clip = s.clip;
            s.source.volume = s.volume;
            s.source.pitch = s.pitch;
            s.source.loop = s.loop;
        }
        
        
        volumeSlider.value = GameSettings.volumeLevel;

        r_Controller.SetupRangeAttackController();

        sceneTransition.SetTrigger("Transition");

        if(currentSceneNumber == 13 && !isKunaiModeTutFinished)
        {
            StartModeTutorial();
        }
       
    }
    

    private void Update()
    {
        if (isGameStart)
        {
            
            
            if (Input.GetButtonDown("Menu") && !r_Controller.isBeganThrowing)
            {
                if (!isGamePaused)
                    OnPause();
                else
                    OnResume();
            }
            if (Input.GetButtonDown("Item"))
            {
                Itemaction();
            }
           
        }
    }


    IEnumerator BeginUITutorial()
    {
        yield return new WaitForSeconds(2);
        Cursor.visible = true;        
        GameObject.FindGameObjectWithTag("Player").GetComponent<Animator>().updateMode = AnimatorUpdateMode.Normal;
        Time.timeScale = 0;
        isGamePaused = true;
        GameSettings.isGamePaused = true;
        TutorialUIElements.SetActive(true);
    }

    public void EndTutorialUI()
    {
        Time.timeScale = 1;
        GameObject.FindGameObjectWithTag("Player").GetComponent<Animator>().updateMode = AnimatorUpdateMode.UnscaledTime;
        Cursor.visible = false;
        isGamePaused = false;
        GameSettings.isGamePaused = false;
    }


    public void BeginKunaiTutorial()
    {
        Cursor.visible = true;
        GameObject.FindGameObjectWithTag("Player").GetComponent<Animator>().updateMode = AnimatorUpdateMode.Normal;
        Time.timeScale = 0;
        isGamePaused = true;
        GameSettings.isGamePaused = true;
        KunaiCoolDownTutorial.SetActive(true);
        tutorialBegin = false;        
    }

    public void GoToNextArea()
    {       

        DataTransferBetweenAreas.unlockedAbilities = new List<int>();
        DataTransferBetweenAreas.equippedAbilitesMap = new List<int>();

        List<int> unlockedAbilitites = new List<int>();
        AbilityManager abilityManager = GetComponent<AbilityManager>();

        for (int i = 0; i < abilityManager.UnlockedAbilities.Length; i++)
        {
            if (abilityManager.UnlockedAbilities[i].GetComponent<Image>().enabled)
            {
                unlockedAbilitites.Add(i);
            }
            else
            {
                unlockedAbilitites.Add(-1);
            }
        }
        DataTransferBetweenAreas.UpdateUnlockedAbilities(unlockedAbilitites);


        List<int> map = new List<int>();
        for (int i = 0; i < abilityManager.EquippedAbilities.Length; i++)
        {

            if (abilityManager.EquippedAbilities[i].GetComponent<Image>().enabled)
            {
                int mapIndex = abilityManager.EquippedAbilitiesDictionary[abilityManager.EquippedAbilities[i]];
                map.Add(mapIndex);
            }

        }

        DataTransferBetweenAreas.UpdateEquippedAbilities(map);


        int sceneToGoInArea = ((4 * currentArea) + 1) - currentSceneNumber;
        int sceneToLoad = currentSceneNumber + sceneToGoInArea;

        StartCoroutine(LoadNextArea(sceneToLoad));
        //SceneManager.LoadScene(5);
    }

    
    IEnumerator LoadNextArea(int scene)
    {
        AsyncOperation operation = SceneManager.LoadSceneAsync(scene);
        while (!operation.isDone)
        {
            yield return null;
        }

        GameSettings.currentScene = scene;
        StartCoroutine(InitiateLevel());
    }
    public void GoToNextRegion()
    {
        currentSceneNumber += 1;

        GameSettings.currentScene = currentSceneNumber;

        StartCoroutine(LoadNextRegion());


       
    }

    IEnumerator LoadNextRegion()
    {
        if((GameSettings.currentScene - 1) % 4 != 0)
            SceneManager.UnloadSceneAsync(GameSettings.currentScene - 1);

        AsyncOperation operation = SceneManager.LoadSceneAsync(GameSettings.currentScene, LoadSceneMode.Additive);
        while (!operation.isDone)
        {
            yield return null;
        }
        string sceneName = SceneManager.GetSceneByBuildIndex(currentSceneNumber).name;
        playerController.transform.position = GameObject.Find(sceneName + "SpawnPoint").transform.position;
        setCameraBounds(currentSceneNumber);
        
        
    }

    #region Methods for Player Health and Stamina
    public void ReducePlayerHealth(float damage)
    {
        float value;
        if (PlayerStats.PlayerHealth >= damage)
        {
            value = PlayerStats.PlayerHealth - damage;
           
        }
        else
        {
            value = 0;
        }
        if (value < 30)
        {
            healthCircleBar.color = Color.red;
        }
        PlayerStats.UpdatePlayerHealth(value);
        //healthBarScript.
        healthBarScript.UpdateHealthBar();

        if (PlayerStats.PlayerHealth <= 0)
        {
            PlayerDied();           
        }
    }
    public void ReducePlayerStamina(float staminaCost)
    {        
        float value = PlayerStats.PlayerStamina - staminaCost;

        PlayerStats.UpdatePlayerStamina(value);

        staminaBarSlider.value = value;
        StartCoroutine(updateStaminaBarOverTime(value));        
    }
    public void IncreasePlayerhealth(float heal)
    {
        float healthMissing = 100 - PlayerStats.PlayerHealth;
        float value = 0;
        if (healthMissing >= heal)
        {
             value = PlayerStats.PlayerHealth + heal;
            
        }
        else
        {
            value = PlayerStats.PlayerHealth + healthMissing;
        }
        if (value > 30)
        {
            healthCircleBar.color = Color.green;
        }
        //StartCoroutine(updateHealthBarOverTime(value));
        healthBarScript.UpdateHealthBar();
        PlayerStats.UpdatePlayerHealth(value);
        //healthBarSlider.value = PlayerStats.PlayerHealth;

    }
    public void IncreasePlayerStamina(float staminaRegen)
    {
        float value;
        if (PlayerStats.PlayerStamina + staminaRegen < 100)
        {
            value = PlayerStats.PlayerStamina + staminaRegen;
           
        }
        else
        {
            value = 100;
        }
        PlayerStats.UpdatePlayerStamina(value);
        staminaBarSlider.value = value;
        StartCoroutine(updateStaminaBarOverTime(value));        
    }
    /*IEnumerator updateHealthBarOverTime(float value)
    {
        while (PlayerStats.PlayerHealth != value)
        {
            PlayerStats.PlayerHealth = Mathf.MoveTowards(PlayerStats.PlayerHealth, value, Time.deltaTime / Time.timeScale);
            yield return null;           
        }
    }*/
    IEnumerator updateStaminaBarOverTime(float value)
    {        
        while (staminaBarWhiteSlider.value != value)
        {
            staminaBarWhiteSlider.value = Mathf.MoveTowards(staminaBarWhiteSlider.value, value, 75 * Time.deltaTime / Time.timeScale);
            yield return null;                     
        }        
    }   
    public void IncreaseMaxPlayerHealth(float value)
    {
        float maxHealth = PlayerStats.PlayerMaxHealth + value;
        PlayerStats.UpdateMaxPlayerHealth(maxHealth);
        healthBarSlider.maxValue = PlayerStats.PlayerMaxHealth;
        healthBarSlider.value = PlayerStats.PlayerHealth;
    }
    public void IncreaseMaxPlayerStamina(float value)
    {
        float maxStealth = PlayerStats.PlayerMaxStamina + value;
        PlayerStats.UpdateMaxPlayerStamina(maxStealth);
        staminaBarSlider.maxValue = PlayerStats.PlayerMaxStamina;
        staminaBarSlider.value = PlayerStats.PlayerStamina;
    }

    public void DecreaseMaxPlayerHealth(float value)
    {
        float maxHealth = PlayerStats.PlayerMaxHealth - value;
        PlayerStats.UpdateMaxPlayerHealth(maxHealth);
        healthBarSlider.maxValue = PlayerStats.PlayerMaxHealth;
        healthBarSlider.value = PlayerStats.PlayerHealth;
    }
    public void DecreaseMaxPlayerStamina(float value)
    {
        float maxStealth = PlayerStats.PlayerMaxStamina - value;
        PlayerStats.UpdateMaxPlayerStamina(maxStealth);
        staminaBarSlider.maxValue = PlayerStats.PlayerMaxStamina;
        staminaBarSlider.value = PlayerStats.PlayerStamina;
    }
   

    #endregion

    #region Getter Methods
    public void NextSceneNumber(int sceneNumber)
    {
        nextSceneNumber = sceneNumber;
    }
    public int GetNextSceneNumber()
    {
        return nextSceneNumber;
    }
    public int GetPreviousSceneNumber()
    {
        return previousSceneNumber;
    }
    public int GetSceneNumber()
    {
        return currentSceneNumber;
    }
    public string GetSceneName()
    {
        return SceneManager.GetSceneByBuildIndex(currentSceneNumber).name;
    }
    #endregion

    

    #region Camera Setup for each Region
    public void setCameraBounds(int sceneNumber)
    {
       
        previousSceneNumber = currentSceneNumber;

        currentSceneNumber = sceneNumber;
        
        CameraController camController = camObject.GetComponent<CameraController>();
        GameObject ob = GameObject.Find(GetSceneName());
        ScenePartLoader loader = ob.GetComponent<ScenePartLoader>();
        Collider2D closestCollider;
        
        if (loader.regions.Length > 1)
        {
            Collider2D[] colliders = new Collider2D[2];
            colliders[0] = loader.regions[0].GetComponent<Collider2D>();
            colliders[1] = loader.regions[1].GetComponent<Collider2D>();


            closestCollider = colliders[0];

            Vector3 closestPointB = closestCollider.ClosestPoint(playerController.transform.position);
            float distanceB = Vector3.Distance(closestPointB, playerController.transform.position);

            foreach (Collider2D collider in colliders)
            {
                Vector3 closestPointA = collider.ClosestPoint(playerController.transform.position);
                float distanceA = Vector3.Distance(closestPointA, playerController.transform.position);

                if (distanceA < distanceB)
                {
                    closestCollider = collider;
                    distanceB = distanceA;
                }
            }
        }
        else
        {
            closestCollider = loader.regions[0].GetComponent<Collider2D>();
            
        }

        camController.border = closestCollider.gameObject;
        int closestPart = 0;
        if(closestCollider.name == loader.regions[0].name)
        {
            closestPart = 0;
        }
        else
        {
            closestPart = 1;
        }

        regionPart = closestPart;

        camController.updateBorder();

        Scene s = SceneManager.GetSceneByBuildIndex(previousSceneNumber);

        GameObject[] gameObjects = s.GetRootGameObjects();
        foreach(GameObject obj in gameObjects)
        {
            if(obj.CompareTag("PathFinder"))
            {
                obj.SetActive(false);
            }
        }
        
        s = SceneManager.GetSceneByBuildIndex(currentSceneNumber);
        gameObjects = s.GetRootGameObjects();
        foreach (GameObject obj in gameObjects)
        {
            if (obj.CompareTag("PathFinder"))
            {
                obj.SetActive(true);
            }


            if (currentSceneNumber > 1)
            {
                if(obj.name == "EnemySpawner")
                {
                    obj.GetComponent<EnemyManager>().enabled = true;
                }               
            }
        }
        GameObject.FindGameObjectWithTag("Player").GetComponent<RangedAttackWithTeleportation>().resetThrowingVars();


    }

    public void switchCameraBoundsBetweenRegion(int part)
    {
        CameraController camController = camObject.GetComponent<CameraController>();
        ScenePartLoader loader = GameObject.Find(GetSceneName()).GetComponent<ScenePartLoader>();

        camController.border = loader.regions[part];
       regionPart = part;
        camController.updateBorder();
    }
    #endregion


    #region Save and Load

    public void DisplaySaves()
    {
        List<string> saveList = SaveSystem.getFileNames();

        if (saveList.Count > 0)
        {
            for (int i = 0; i < saveList.Count; i++)
            {
                saveslotsTexts[i].text = saveList[i];
            }
        }
    }

    public void SaveGame(int slot)
    {
        PlayerData data = new PlayerData();
        data.UpdatePlayerPosition(playerController.transform.position);


        int currentArea = (currentSceneNumber - 1) / 4 + 1;

        data.UpdateArea(currentArea);
        data.updateRegion(currentSceneNumber);


        /*        
        List<int> items = new List<int>();
        Inventory inventory = transform.GetComponent<Inventory>();

        for (int i = 0; i < inventory.slots.Length; i++)
        {
            if (inventory.isFull[i])
            {
                GameObject child = inventory.slots[i].transform.GetChild(0).gameObject;
                switch (child.GetComponent<PickUp>().itemType)
                {
                    case PickUp.ItemType.HealthPotion:
                        items.Add(0);
                        break;
                    case PickUp.ItemType.StaminaPotion:
                        //items.Add(1);
                        break;
                    default:
                        break;
                }        
            }

        }
        data.updateInventory(items);
        
         */

        List<int> unlockedAbilitites = new List<int>();
        AbilityManager abilityManager = GetComponent<AbilityManager>();

        for(int i = 0; i < abilityManager.UnlockedAbilities.Length; i++)
        {
            if(abilityManager.UnlockedAbilities[i].GetComponent<Image>().enabled)
            {
                unlockedAbilitites.Add(i);
            }
            else
            {
                unlockedAbilitites.Add(-1);
            }
        }
        data.UpdateUnlockedAbilities(unlockedAbilitites);

        
        List<int> map = new List<int>();
        for (int i = 0; i < abilityManager.EquippedAbilities.Length; i++)
        {

            if (abilityManager.EquippedAbilities[i].GetComponent<Image>().enabled )
            {
                int mapIndex = abilityManager.EquippedAbilitiesDictionary[abilityManager.EquippedAbilities[i]];
                map.Add(mapIndex);
            }
          
        }

        data.UpdateEquippedAbilities(map);

        data.UpdateCoinCount(coinsCollected);
        
        SaveSystem.SaveGame(data, slot);
        SaveData saveData = new SaveData(data);
        GameSettings.saves.Add(saveData);
        if (isGamePaused)
            OnResume();
        GameSettings.saveSlot = slot;
        GameSettings.currentArea = currentArea;
        GameSettings.currentScene = currentSceneNumber;
    }

    public void LoadGame(int slot)
    {
        SaveData data = SaveSystem.LoadGame(slot);
        currentSceneNumber = data.currentRegion;

        PlayerStats.UpdatePlayerHealth(data.playerHealth);
        PlayerStats.UpdatePlayerStamina(data.playerStamina);

        IncreasePlayerhealth(0);
        IncreasePlayerStamina(0);
        
        setCameraBounds(data.currentRegion);

        Vector3 pos = new Vector3(data.playerPosition[0], data.playerPosition[1], data.playerPosition[2]);

        playerController.transform.position = pos;     

        AbilityManager manager = GetComponent<AbilityManager>();

        manager.Setup();

        for(int i =0; i <data.unlockedAbilities.Count; i++)
        {
            if (data.unlockedAbilities[i] != -1)
            {
                int index = data.unlockedAbilities[i];
                GameObject obj = manager.UnlockedAbilities[index];
                manager.ReUnlockAbility(obj.name);
            }
        }

        for(int i = 0; i < data.equippedAbilititesMap.Count; i++)
        {
            
            int mapIndex = data.equippedAbilititesMap[i];
            if (mapIndex != -1)
            {
                manager.equipAbility(mapIndex);
            }
        }


        UpdateCoin(data.coinCount);

        GameSettings.currentArea = data.currentArea;
        GameSettings.currentScene = data.currentRegion;
        GameSettings.playerPos = data.playerPosition;

        Sound[] sounds = FindObjectOfType<AudioManager>().sounds;
        foreach (Sound s in sounds)
        {
            s.source = gameObject.AddComponent<AudioSource>();
            s.source.clip = s.clip;
            s.source.volume = s.volume;
            s.source.pitch = s.pitch;
            s.source.loop = s.loop;
        }
        volumeSlider.value = data.gameVolumeLevel;
        VolumeChanged();

        if(currentSceneNumber == 13)
        {
            isKunaiModeTutFinished = true;
        }

    }



    #endregion



    public int GetRegionPart()
    {
        return regionPart;
    }

    // Player Died, Reload
    public void PlayerDied()
    {
        playerController.PlayerDied();
    }

  
    public void Itemaction()
    {
        if (isInventoryOpen == false)
        {            
            Inventory.SetActive(true);
            isInventoryOpen = true;
            Cursor.visible = true;
            Time.timeScale = 0;
            isGamePaused = true;
            GameSettings.isGamePaused = true;
        }
        else
        {
            
            Inventory.SetActive(false);
            isInventoryOpen = false;
            Cursor.visible = false;
            Time.timeScale = 1f;
            isGamePaused = false;
            GameSettings.isGamePaused = false;
        }

    }

    public void InstantiateHealthPotion(Vector3 position)
    {
        int random = Random.Range(0, 11);
        if (random < 5)
        {
            GameObject potion =  Instantiate(healthPotion, position, Quaternion.identity);
            potion.GetComponent<Rigidbody2D>().AddForce(Vector2.up * 50);
        }
    }



    public void EnableAbilityInGame(string abilityName)
    {
        switch (abilityName)
        {
            case "DashAbility":
                playerController.dashUnlocked = true;
                break;

            case "SpecialDash":
                playerController.dashUnlocked = true;
                playerController.dashUpgraded = true;
                break;

            case "WallThrow":
                playerController.wallThrowUnlocked = true;
                break;
            case "SpellAbility":
                playerController.spellAbilityUnlocked = true;
                break;
                    


            default:
                break;
        }
        
    }
    
    public void DisableAbilityInGame(string abilityName)
    {
        switch (abilityName)
        {
            case "DashAbility":
                playerController.dashUnlocked = false;
                break;

            case "SpecialDash":
                playerController.dashUnlocked = false;
                playerController.dashUpgraded = false;
                break;

            case "WallThrow":
                playerController.wallThrowUnlocked = false;
                break;
            case "SpellAbility":
                playerController.spellAbilityUnlocked = false;
                break;

            default:
                break;
        }
    }


    public void UnlockAbility(string abilityName)
    {
        AbilityManager abilityManager = GetComponent<AbilityManager>();
        abilityManager.UnlockAbility(abilityName);
    }

    public void AddCoin()
    {
        if (coinsCollected < maxCoinsAllowed)
        {
            coinsCollected++;
            coinCount.text = coinsCollected.ToString();            
        }

    }

    public void UpdateCoin(int count)
    {
        coinsCollected = count;
        coinCount.text = coinsCollected.ToString();
    }


    public void OnPause()
    {
        Cursor.visible = true;
        Time.timeScale = 0;
        isGamePaused = true;
        GameSettings.isGamePaused = true;
        ingameMenu.SetActive(true);
        upgradeManager.enabled = true;

    }


    public void OnResume()
    {
        Cursor.visible = false;
        Time.timeScale = 1f;
        isGamePaused = false;
        GameSettings.isGamePaused = false;
        ingameMenu.SetActive(false);
        upgradeManager.enabled = false;
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

    public void OpenSaveMenu()
    {
        
        Cursor.visible = true;
        Time.timeScale = 0;
        isGamePaused = true;
        GameSettings.isGamePaused = true;
        saveMenu.SetActive(true);
        DisplaySaves();
        upgradeManager.enabled = true;
    }

    public void ReturnToMainMenu()
    {
        SceneManager.LoadScene(0);
    }


    public void FinishKunaiTutorial()
    {
        r_Controller.canThrow = true;
    }

    private void StartModeTutorial()
    {
        Cursor.visible = true;
        //GameObject.FindGameObjectWithTag("Player").GetComponent<Animator>().updateMode = AnimatorUpdateMode.Normal;
        //Time.timeScale = 0;
        isGamePaused = true;
        GameSettings.isGamePaused = true;
        mode2.SetActive(true);
        kunaiModeTutObj.SetActive(true);
    }
    public void EnableKunaiModes()
    {       
        isKunaiModeUnlocked = true;
        isGamePaused = false;
        GameSettings.isGamePaused = false;
    }
}
