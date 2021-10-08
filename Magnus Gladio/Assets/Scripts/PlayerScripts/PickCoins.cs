using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class PickCoins : MonoBehaviour
{
    private void OnTriggerEnter2D(Collider2D collision)
    {
        if(collision.CompareTag("PlayerCircleCollider"))
        {
            GameObject.FindGameObjectWithTag("Game Manager").GetComponent<GameManager>().AddCoin();
            Destroy(gameObject);
        }
    }
}
