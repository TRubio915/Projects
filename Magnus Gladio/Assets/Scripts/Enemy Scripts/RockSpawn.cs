using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class RockSpawn : MonoBehaviour
{
    public GameObject rockPrefab; //Prefab for rock spawn
    public float spawnTime;
    public float lifetime; //Lifetime of rock prefab
    // Start is called before the first frame update
    void Start()
    {
        //Invoke spawn after lifetime, and every lifetime x2
        InvokeRepeating("SpawnRock", spawnTime, (spawnTime * 2));
    }

    private void SpawnRock()
    {
        GameObject r = Instantiate(rockPrefab, transform.position, Quaternion.identity);
        Destroy(r, lifetime);
    }
}
