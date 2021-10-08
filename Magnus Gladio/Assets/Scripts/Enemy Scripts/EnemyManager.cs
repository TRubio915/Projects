using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.SceneManagement;

public class EnemyManager : MonoBehaviour
{

    public int region , area;

    private GameObject spawner;
    private Transform[] enemySpawnPoints;

    private bool enemiesSpawned = false;
    
    [SerializeField]
    private List<GameObject> enemyPrefabs;
    
    private GameObject[] enemies;

    // Start is called before the first frame update
    void Start()
    {
        if (!enemiesSpawned)
        {
            foreach (GameObject obj in SceneManager.GetSceneByName("Area"+area+"Region" + region).GetRootGameObjects())
            {
                if (obj.name == "EnemySpawnPoints")
                {
                    spawner = obj;
                    break;
                }
            }

            enemySpawnPoints = spawner.transform.GetComponentsInChildren<Transform>();
            enemies = new GameObject[enemySpawnPoints.Length];
            for (int i = 0; i < enemySpawnPoints.Length; i++)
            {
                int randomIndex = Random.Range(0, 2);
                enemies[i] = Instantiate(enemyPrefabs[randomIndex], enemySpawnPoints[i].position, Quaternion.identity);
                //enemies[i].transform.parent = this.transform;
            }
            enemiesSpawned = true;
        }
    }

    public void UnloadEnemies()
    {
        if (enemiesSpawned)
        {
            for (int i = 0; i < enemies.Length; i++)
            {
                Destroy(enemies[i]);
            }
        }
    }
}
