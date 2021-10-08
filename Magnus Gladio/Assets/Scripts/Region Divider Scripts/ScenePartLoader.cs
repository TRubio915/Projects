using UnityEngine;
using UnityEngine.SceneManagement;

public enum CheckMethod
{
    Distance,
    Trigger
}

public class ScenePartLoader : MonoBehaviour
{
    public int sceneNumber;
    public Transform player;
    public CheckMethod checkMethod;
    public float loadRange;

    //Scene state
    private bool isLoaded;
    private bool shouldLoad;
    
    public GameObject[] regions;

    [SerializeField]
    private GameManager manager;

    
   

    void Start()
    {
        manager = GameObject.FindGameObjectWithTag("Game Manager").GetComponent<GameManager>();
        player = GameObject.FindGameObjectWithTag("Player").transform;
        //verify if the scene is already open to avoid opening a scene twice
        if (SceneManager.sceneCount > 0)
        {
            for (int i = 0; i < SceneManager.sceneCount; ++i)
            {
                Scene scene = SceneManager.GetSceneAt(i);
                if (scene.name == gameObject.name)
                {
                    isLoaded = true;
                }
            }
        }
    }

    void Update()
    {
        //Checking which method to use
        if (checkMethod == CheckMethod.Distance)
        {
            DistanceCheck();
        }
        else if (checkMethod == CheckMethod.Trigger)
        {
            TriggerCheck();
        }
    }

    void DistanceCheck()
    {
        //Checking if the player is within the range
        if (Vector3.Distance(player.position, transform.position) < loadRange)
        {
            LoadScene();
        }
        else
        {
            UnLoadScene();
        }
    }

    void LoadScene()
    {
        if (!isLoaded && !SceneManager.GetSceneByName(gameObject.name).isLoaded)
        {
            //Loading the scene, using the gameobject name as it's the same as the name of the scene to load
            SceneManager.LoadSceneAsync(gameObject.name, LoadSceneMode.Additive);
            //We set it to true to avoid loading the scene twice
            isLoaded = true;
            sceneNumber = SceneManager.GetSceneByName(gameObject.name).buildIndex;
            manager.NextSceneNumber(sceneNumber);
            //manager.setCameraBounds();
        }
    }

    void UnLoadScene()
    {
        if (isLoaded)
        {
            foreach (GameObject obj in SceneManager.GetSceneByName(gameObject.name).GetRootGameObjects())
            {
                if (obj.name == "EnemySpawner")
                {
                    obj.GetComponent<EnemyManager>().UnloadEnemies();
                }
            }
            
            SceneManager.UnloadSceneAsync(gameObject.name);
            isLoaded = false;
        }
    }

    private void OnTriggerEnter2D(Collider2D other)
    {
        if (other.CompareTag("Player"))
        {
            shouldLoad = true;
        }
    }

    private void OnTriggerExit2D(Collider2D other)
    {
        if (other.CompareTag("Player"))
        {
            shouldLoad = false;
        }
    }

    void TriggerCheck()
    {
        //shouldLoad is set from the Trigger methods
        if (shouldLoad)
        {
            LoadScene();
        }
        else
        {
            UnLoadScene();
        }
    }
    private void OnDrawGizmos()
    {
        Gizmos.DrawWireSphere(transform.position, loadRange);
        //Gizmos.DrawWireCube(transform.position, GetComponent<PolygonCollider2D>().bounds.size);
    }




}